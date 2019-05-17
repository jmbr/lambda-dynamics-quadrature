cmake_minimum_required(VERSION 3.1)

project(lambda-dynamics-quadrature C CXX)
enable_language(C CXX)

find_package(Boost REQUIRED)
find_package(Catch2)

set(CXX_WARN_FLAGS "-W -Wall -Wextra -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")
set(CMAKE_CXX_FLAGS "-std=c++17 ${CXX_WARN_FLAGS}")
set(CMAKE_MACOSX_RPATH 1)

find_program(RUN_CLANG_TIDY run-clang-tidy.py DOC "Python script to run clang-tidy")
if(RUN_CLANG_TIDY)
  file(GLOB_RECURSE SOURCE_FILES *.cpp *.h)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  add_custom_target(lint
    COMMAND run-clang-tidy.py -header-filter='${CMAKE_SOURCE_DIR}/.*'
                              -p=${CMAKE_BINARY_DIR} ${SOURCE_FILES})
endif()

add_subdirectory(contrib/cubature EXCLUDE_FROM_ALL)
include_directories(${CMAKE_SOURCE_DIR}/contrib/cubature ${Boost_INCLUDE_DIRS})
link_directories(${cubature_BINARY_DIR})

add_library(lambda-dynamics-quadrature potential.cpp ensemble-average.cpp)

add_executable(run-computation run-computation.cpp)
target_link_libraries(run-computation lambda-dynamics-quadrature cubature ${Boost_LIBRARIES})

if(Catch2_FOUND)
  add_executable(test-main
    test-main.cpp
    test-ensemble-average.cpp
    test-potential.cpp
    test-quadrature.cpp)

  target_link_libraries(test-main lambda-dynamics-quadrature cubature ${Boost_LIBRARIES} Catch2::Catch2)

  include(CTest)
  include(Catch)
  catch_discover_tests(test-main)
endif()
