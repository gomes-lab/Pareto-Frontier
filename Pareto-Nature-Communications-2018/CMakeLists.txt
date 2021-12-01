cmake_minimum_required(VERSION 3.3.2)
project(Amazon-E-GHG)

set(CMAKE_CXX_STANDARD 11)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(SOURCE_FILES
        DP_Algorithm.cpp
        DP_Algorithm.h
        HyperNet.cpp
        HyperNet.h
        main.cpp
        Pareto_Solution.cpp
        Pareto_Solution.h
        ThreadPool.h)

add_executable(Amazon-E-GHG ${SOURCE_FILES})
target_link_libraries(Amazon-E-GHG Threads::Threads)
