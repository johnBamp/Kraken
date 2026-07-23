cmake_minimum_required(VERSION 3.23)

project(Kraken LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_library(Kraken INTERFACE)

target_include_directories(Kraken
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include

)

add_subdirectory(examples)