cmake_minimum_required(VERSION 3.12)
project(WHFC)

find_package(TBB REQUIRED)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wunused")
set(CMAKE_CXX_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lm")
#set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer")

add_executable(WHFC main.cpp)
target_link_libraries(WHFC PUBLIC TBB::tbb TBB::tbbmalloc)

add_executable(TESTS run_tests.cpp)
target_link_libraries(TESTS PUBLIC TBB::tbb TBB::tbbmalloc)

add_executable(SnapshotTester snapshot_tester.cpp)
target_link_libraries(SnapshotTester PUBLIC TBB::tbb TBB::tbbmalloc)

add_executable(FlowTester flow_tester.cpp)
target_link_libraries(FlowTester PUBLIC TBB::tbb TBB::tbbmalloc)
