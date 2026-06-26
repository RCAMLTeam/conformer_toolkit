"""Automatic K selection for KNN graphs built from conformer RMSD matrices."""

from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np


@dataclass(frozen=True)
class KNNGraph:
    """Undirected KNN graph built from an RMSD matrix."""

    k: int
    edges: np.ndarray
    edge_weights: np.ndarray
    neighbor_indices: np.ndarray

    @property
    def edge_count(self) -> int:
        return int(self.edges.shape[0])

    @property
    def avg_edge_weight(self) -> float:
        if self.edge_weights.size == 0:
            return float("nan")
        return float(np.mean(self.edge_weights))


@dataclass(frozen=True)
class KDiagnostic:
    """Connectivity and quality diagnostics for one tested K value."""

    k: int
    connected: bool
    n_components: int
    edge_count: int
    avg_edge_weight: float
    connectivity_score: float
    stability_score: float


@dataclass(frozen=True)
class KOptimizationResult:
    """Result from smallest-connected-K search."""

    optimal_k: int
    connected: bool
    avg_edge_weight: float
    edge_count: int
    diagnostics: list[KDiagnostic]
    warning: str | None = None


class UnionFind:
    """Small Union-Find implementation for repeated graph connectivity checks."""

    def __init__(self, n_items: int) -> None:
        self.parent = np.arange(n_items, dtype=np.int64)
        self.rank = np.zeros(n_items, dtype=np.int8)
        self.components = n_items

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = int(self.parent[root])
        while self.parent[item] != item:
            parent = int(self.parent[item])
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return

        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1
        self.components -= 1


class KNNKOptimizer:
    """Find the smallest connected K for conformer KNN graphs.

    Parameters
    ----------
    rmsd_matrix:
        Square precomputed RMSD matrix. The diagonal is ignored.
    k_min:
        First K to test. Defaults to 2.
    k_max_cap:
        Maximum search cap before applying the N - 1 bound. Defaults to 50.
    edge_mode:
        ``"one_way"`` adds an undirected edge when either endpoint selects the
        other as a neighbor. ``"mutual"`` only keeps reciprocal neighbor pairs.
    symmetrize:
        If true, use ``0.5 * (matrix + matrix.T)`` for edge weights/selection.
        RDKit RMSD matrices should already be symmetric; this just removes tiny
        numerical asymmetries.
    """

    def __init__(
        self,
        rmsd_matrix: np.ndarray,
        *,
        k_min: int = 2,
        k_max_cap: int = 50,
        edge_mode: str = "one_way",
        symmetrize: bool = True,
    ) -> None:
        matrix = np.asarray(rmsd_matrix, dtype=np.float64)
        self._validate_matrix(matrix)
        if edge_mode not in {"one_way", "mutual"}:
            raise ValueError("edge_mode must be 'one_way' or 'mutual'")
        if k_min < 1:
            raise ValueError("k_min must be at least 1")
        if k_max_cap < 1:
            raise ValueError("k_max_cap must be at least 1")

        self.n_nodes = int(matrix.shape[0])
        self.k_min = int(k_min)
        self.k_max_cap = int(k_max_cap)
        self.edge_mode = edge_mode
        self.distances = 0.5 * (matrix + matrix.T) if symmetrize else matrix.copy()
        np.fill_diagonal(self.distances, np.inf)

    @staticmethod
    def _validate_matrix(matrix: np.ndarray) -> None:
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("rmsd_matrix must be a square N x N array")
        if matrix.shape[0] < 2:
            raise ValueError("rmsd_matrix must contain at least two conformers")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("rmsd_matrix must contain only finite values")
        if np.any(matrix < 0.0):
            raise ValueError("rmsd_matrix cannot contain negative distances")

    @property
    def k_max(self) -> int:
        return min(self.n_nodes - 1, self.k_max_cap)

    def build_graph(self, k: int) -> KNNGraph:
        """Build an undirected KNN graph for K."""

        if k < 1 or k > self.n_nodes - 1:
            raise ValueError(f"k must be in [1, {self.n_nodes - 1}]")

        neighbor_indices = self._knn_indices(k)
        edges = self._edges_from_neighbors(neighbor_indices)
        edge_weights = self.distances[edges[:, 0], edges[:, 1]] if edges.size else np.array([], dtype=np.float64)
        return KNNGraph(k=k, edges=edges, edge_weights=edge_weights, neighbor_indices=neighbor_indices)

    def optimize(self) -> KOptimizationResult:
        """Return the first K in [k_min, min(N - 1, 50)] that is connected."""

        diagnostics: list[KDiagnostic] = []
        last_graph: KNNGraph | None = None
        start_k = min(self.k_min, self.k_max)
        for k in range(start_k, self.k_max + 1):
            graph = self.build_graph(k)
            last_graph = graph
            diagnostic = self.diagnose_graph(graph)
            diagnostics.append(diagnostic)
            if diagnostic.connected:
                return KOptimizationResult(
                    optimal_k=k,
                    connected=True,
                    avg_edge_weight=diagnostic.avg_edge_weight,
                    edge_count=diagnostic.edge_count,
                    diagnostics=diagnostics,
                )

        assert last_graph is not None
        warning = f"no connected KNN graph found for K in [{self.k_min}, {self.k_max}]"
        warnings.warn(warning, RuntimeWarning, stacklevel=2)
        last_diagnostic = diagnostics[-1]
        return KOptimizationResult(
            optimal_k=self.k_max,
            connected=False,
            avg_edge_weight=last_diagnostic.avg_edge_weight,
            edge_count=last_diagnostic.edge_count,
            diagnostics=diagnostics,
            warning=warning,
        )

    def diagnose_graph(self, graph: KNNGraph) -> KDiagnostic:
        n_components = self._component_count(graph.edges)
        connected = n_components == 1
        return KDiagnostic(
            k=graph.k,
            connected=connected,
            n_components=n_components,
            edge_count=graph.edge_count,
            avg_edge_weight=graph.avg_edge_weight,
            connectivity_score=self.connectivity_score(n_components),
            stability_score=self.stability_score(graph.neighbor_indices, graph.edges),
        )

    @staticmethod
    def connectivity_score(n_components: int) -> float:
        """Strategy 4 connectivity term: 1 if connected, else 1 / n_components."""

        if n_components < 1:
            raise ValueError("n_components must be positive")
        return 1.0 / float(n_components)

    @staticmethod
    def stability_score(neighbor_indices: np.ndarray, edges: np.ndarray) -> float:
        """Mean Jaccard similarity of endpoint neighbor sets across graph edges."""

        if edges.size == 0:
            return 0.0
        neighbor_sets = [set(row.tolist()) for row in neighbor_indices]
        similarities = []
        for left, right in edges:
            left_neighbors = neighbor_sets[int(left)]
            right_neighbors = neighbor_sets[int(right)]
            union = left_neighbors | right_neighbors
            similarities.append(len(left_neighbors & right_neighbors) / len(union) if union else 0.0)
        return float(np.mean(similarities))

    @staticmethod
    def strategy4_score(
        diagnostic: KDiagnostic,
        *,
        alpha: float = 1.0,
        beta: float = 0.01,
        gamma: float = 0.1,
    ) -> float:
        """Weighted Strategy 4 score hook for future K selection policies.

        The default signs make higher better: connectivity and stability add to
        the score, while larger average RMSD is penalized. Keep alpha dominant
        when connectivity is the primary constraint.
        """

        return (
            alpha * diagnostic.connectivity_score
            - beta * diagnostic.avg_edge_weight
            + gamma * diagnostic.stability_score
        )

    def _knn_indices(self, k: int) -> np.ndarray:
        candidates = np.argpartition(self.distances, kth=k - 1, axis=1)[:, :k]
        candidate_distances = np.take_along_axis(self.distances, candidates, axis=1)
        order = np.argsort(candidate_distances, axis=1)
        return np.take_along_axis(candidates, order, axis=1)

    def _edges_from_neighbors(self, neighbor_indices: np.ndarray) -> np.ndarray:
        source = np.repeat(np.arange(self.n_nodes, dtype=np.int64), neighbor_indices.shape[1])
        target = neighbor_indices.reshape(-1).astype(np.int64, copy=False)
        pairs = np.column_stack((np.minimum(source, target), np.maximum(source, target)))
        pairs = pairs[pairs[:, 0] != pairs[:, 1]]
        if pairs.size == 0:
            return np.empty((0, 2), dtype=np.int64)

        edges, counts = np.unique(pairs, axis=0, return_counts=True)
        if self.edge_mode == "mutual":
            edges = edges[counts == 2]
        return edges.astype(np.int64, copy=False)

    def _component_count(self, edges: np.ndarray) -> int:
        union_find = UnionFind(self.n_nodes)
        for left, right in edges:
            union_find.union(int(left), int(right))
        return union_find.components
