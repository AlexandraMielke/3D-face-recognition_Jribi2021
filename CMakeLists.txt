/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */
# CMakeList.txt: CMake-Projekt für "2021_Jribi". Schließen Sie die Quelle ein, und definieren Sie
# projektspezifische Logik hier.
#
cmake_minimum_required (VERSION 3.8)
set( CMAKE_CXX_COMPILER "/usr/bin/g++" )
set( CMAKE_C_COMPILER "/usr/bin/gcc" )

project("JribiDescription")

#Source files
file(GLOB_RECURSE lib_sources "lib/SplineLib/*.cpp")
set(sources ${sources} "src/2021Jribi.cpp" "src/LevelCurves.cpp" "src/pywrap.cpp"
                        ${lib_sources})

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(JribiDescription ${sources})

#target definition only needed when pybind is not used
#add_executable(2021_Jribi ${sources})


#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET 2021_Jribi PROPERTY CXX_STANDARD 20)
#endif()

# TODO: Fügen Sie bei Bedarf Tests hinzu, und installieren Sie Ziele.
