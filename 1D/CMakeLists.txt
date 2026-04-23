cmake_minimum_required(VERSION 3.10)

project(Heat1D VERSION 1.0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(EXE "Heat1D_exe")

add_executable(${EXE}
    Initial_data.cc
    Main.cc
)

find_package(HDF5 REQUIRED COMPONENTS C)
target_link_libraries(${EXE} PRIVATE HDF5::HDF5)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

install(TARGETS ${EXE}
        RUNTIME DESTINATION bin
)
