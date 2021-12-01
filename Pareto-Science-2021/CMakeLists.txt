cmake_minimum_required(VERSION 3.3.2)
project(Amazon_DP)

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
        Pareto_Opt_List.cpp
        Pareto_Opt_List.h
        Pareto_Opt_Node.cpp 
        Pareto_Opt_Node.h
        Frontier_List.cpp
        Frontier_List.h
        ThreadPool.h)

add_executable(Amazon ${SOURCE_FILES})
target_link_libraries(Amazon Threads::Threads)