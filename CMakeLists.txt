cmake_minimum_required(VERSION 3.0)
project(hexex2geogram)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(ultimaille)

include_directories(ultimaille ultimaille/ext)

set(LIBS
    ultimaille
)

add_executable(hexex2geogram main.cpp)
target_link_libraries(hexex2geogram ${LIBS})