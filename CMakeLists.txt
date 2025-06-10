cmake_minimum_required(VERSION 3.16)
project(LinearAlgebraSolver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(project
        main.cpp
        Vector.cpp
        Matrix.cpp
        LinearSystem.cpp
        PosSymLinSystem.cpp
        GeneralLinearSystem.cpp
)
