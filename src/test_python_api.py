from pathlib import Path

from conformer_toolkit import Conformer_Group


PATH = Path(__file__).parent / "testdata" / "h2_energies.xyz"
TEMPLATE = "frame = {frame} energy = {energy}"


def load():
    return Conformer_Group.from_multi_xyz(str(PATH), comment_template=TEMPLATE)


group = load()
assert group.records()[0].properties == {"frame": "0", "energy": "3.0"}

placeholder_group = Conformer_Group.from_multi_xyz(
    str(PATH), comment_template="frame = {_framenumber} energy = {energy}"
)
assert placeholder_group.records()[0].properties == {"energy": "3.0"}

group.sort_by_energy()
assert [record.properties["frame"] for record in group.records()] == ["1", "3", "2", "0"]

group = load()
group.filter_by_maximum_energy(1.5)
assert [record.properties["frame"] for record in group.records()] == ["1", "3"]

group = load()
group.retain_lowest_energy_percent(50)
assert [record.properties["frame"] for record in group.records()] == ["1", "3"]

group = load()
group.filter_by_boltzmann_population_ratio(0.5)
assert [record.properties["frame"] for record in group.records()] == ["1", "3"]
