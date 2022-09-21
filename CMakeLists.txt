cmake_minimum_required(VERSION 3.0)
project(hexex2geogram)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if (NOT WIN32) #For ultimaille
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  add_compile_definitions(USE_OPENMP true)
endif()
endif()

add_subdirectory(ultimaille)

include_directories(ultimaille ultimaille/ext)

set(LIBS
    ultimaille
)

add_executable(hexex2geogram main.cpp)
target_link_libraries(hexex2geogram ${LIBS})