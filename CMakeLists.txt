cmake_minimum_required(VERSION 3.12)
project(MorseCodeMines VERSION 1.0)

# Specify C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable debug symbols
set(CMAKE_BUILD_TYPE Debug)

# Add compiler-specific flags (equivalent to -g)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-g)
endif()

# Add all source files in the project
file(GLOB_RECURSE SOURCES 
    "*.cpp"
)
# Create an executable for each cpp file
foreach(SOURCE_FILE ${SOURCES})
    get_filename_component(FILE_NAME ${SOURCE_FILE} NAME_WE)
    add_executable(${FILE_NAME} ${SOURCE_FILE})
endforeach()

# For specific executable (if you want to build just the current file)
# add_executable(decoding_v1 decoding_v1.cpp)