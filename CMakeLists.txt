cmake_minimum_required(VERSION 3.18)
project(conformer_toolkit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(RDKIT_ROOT "" CACHE PATH "RDKit installation prefix")
set(RDKIT_SOURCE_ROOT "" CACHE PATH "RDKit source tree, needed when DetermineBonds is not installed")

set(_rdkit_include_hints
    "${RDKIT_ROOT}/include/rdkit"
    "${RDKIT_ROOT}/usr/include/rdkit"
    "/usr/include/rdkit"
    "/usr/local/include/rdkit"
)
set(_rdkit_library_hints
    "${RDKIT_ROOT}/lib"
    "${RDKIT_ROOT}/usr/lib"
    "${RDKIT_ROOT}/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    "${RDKIT_ROOT}/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
)

find_path(RDKIT_INCLUDE_DIR RDGeneral/export.h PATHS ${_rdkit_include_hints})
if(NOT RDKIT_INCLUDE_DIR)
    message(FATAL_ERROR "RDKit headers not found. Set RDKIT_ROOT.")
endif()

function(find_rdkit_library output name)
    find_library(${output} NAMES "RDKit${name}" PATHS ${_rdkit_library_hints})
    if(NOT ${output})
        message(FATAL_ERROR "RDKit library RDKit${name} not found. Set RDKIT_ROOT.")
    endif()
    set(${output} "${${output}}" PARENT_SCOPE)
endfunction()

find_rdkit_library(RDKIT_FILEPARSERS FileParsers)
find_rdkit_library(RDKIT_SUBSTRUCTMATCH SubstructMatch)
find_rdkit_library(RDKIT_GRAPHMOL GraphMol)
find_rdkit_library(RDKIT_MOLALIGN MolAlign)
find_rdkit_library(RDKIT_GEOMETRY RDGeometryLib)
find_rdkit_library(RDKIT_GENERAL RDGeneral)

set(RDKIT_LIBRARIES
    ${RDKIT_FILEPARSERS}
    ${RDKIT_SUBSTRUCTMATCH}
    ${RDKIT_GRAPHMOL}
    ${RDKIT_MOLALIGN}
    ${RDKIT_GEOMETRY}
    ${RDKIT_GENERAL}
)

find_path(
    RDKIT_DETERMINEBONDS_INCLUDE_DIR
    GraphMol/DetermineBonds/DetermineBonds.h
    PATHS ${_rdkit_include_hints}
)
find_library(
    RDKIT_DETERMINEBONDS_LIBRARY
    NAMES RDKitDetermineBonds
    PATHS ${_rdkit_library_hints}
)

if(RDKIT_DETERMINEBONDS_INCLUDE_DIR AND RDKIT_DETERMINEBONDS_LIBRARY)
    list(APPEND RDKIT_LIBRARIES ${RDKIT_DETERMINEBONDS_LIBRARY})
    set(DETERMINEBONDS_SOURCES)
    set(DETERMINEBONDS_INCLUDE_DIR "${RDKIT_DETERMINEBONDS_INCLUDE_DIR}")
else()
    if(NOT EXISTS "${RDKIT_SOURCE_ROOT}/Code/GraphMol/DetermineBonds/DetermineBonds.cpp")
        message(FATAL_ERROR
            "This RDKit installation omits DetermineBonds. "
            "Set RDKIT_SOURCE_ROOT to a matching RDKit source checkout."
        )
    endif()
    set(DETERMINEBONDS_SOURCES
        "${RDKIT_SOURCE_ROOT}/Code/GraphMol/DetermineBonds/DetermineBonds.cpp"
    )
    set(DETERMINEBONDS_INCLUDE_DIR "${RDKIT_SOURCE_ROOT}/Code")
endif()

function(configure_conformer_target target)
    target_include_directories(${target} PRIVATE
        "${RDKIT_INCLUDE_DIR}"
        "${DETERMINEBONDS_INCLUDE_DIR}"
        "${RDKIT_ROOT}/include"
        "${RDKIT_ROOT}/usr/include"
    )
    target_link_libraries(${target} PRIVATE ${RDKIT_LIBRARIES})
    if(UNIX AND NOT APPLE)
        # DT_RPATH is inherited by transitive RDKit dependencies. DT_RUNPATH is
        # not, which breaks repository-local SDKs whose libraries are not also
        # installed in the system loader path.
        target_link_options(${target} PRIVATE "-Wl,--disable-new-dtags")
    endif()
    target_compile_definitions(${target} PRIVATE RDKIT_DETERMINEBONDS_BUILD)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    set_target_properties(${target} PROPERTIES
        BUILD_RPATH "${_rdkit_library_hints}"
        INSTALL_RPATH "${_rdkit_library_hints}"
    )
endfunction()

add_executable(
    conformer_deduplicate
    src/conformer_deduplicate.cpp
    ${DETERMINEBONDS_SOURCES}
)
configure_conformer_target(conformer_deduplicate)

add_executable(
    conformer_identical
    src/conformer_identical.cpp
    ${DETERMINEBONDS_SOURCES}
)
configure_conformer_target(conformer_identical)

option(BUILD_PYTHON_BINDINGS "Build the pybind11 extension" ON)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Development.Module REQUIRED)
    find_path(
        PYBIND11_INCLUDE_DIR
        pybind11/pybind11.h
        PATHS
            "${RDKIT_ROOT}/include"
            "${RDKIT_ROOT}/usr/include"
            "/usr/include"
            "/usr/local/include"
    )
    if(NOT PYBIND11_INCLUDE_DIR)
        message(FATAL_ERROR "pybind11 headers not found. Install pybind11-dev or set RDKIT_ROOT to the local SDK.")
    endif()

    Python3_add_library(
        conformer_toolkit_cpp
        MODULE
        WITH_SOABI
        src/conformer_toolkit_bindings.cpp
        ${DETERMINEBONDS_SOURCES}
    )
    target_include_directories(conformer_toolkit_cpp PRIVATE "${PYBIND11_INCLUDE_DIR}")
    configure_conformer_target(conformer_toolkit_cpp)
endif()

set_target_properties(conformer_deduplicate conformer_identical PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src"
)
if(TARGET conformer_toolkit_cpp)
    set_target_properties(conformer_toolkit_cpp PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src"
    )
endif()

enable_testing()
add_test(
    NAME conformer_regressions
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/src/test_cpp.py"
        --identical "$<TARGET_FILE:conformer_identical>"
        --deduplicate "$<TARGET_FILE:conformer_deduplicate>"
)
