cmake_minimum_required(VERSION 3.24)
project(myMPM)

set(CMAKE_CXX_STANDARD 17)

# 头文件路径，用在 add_executable() 之前
include_directories(external/partio/src/lib)
include_directories(external/spdlog/include)
include_directories(include)

add_executable(myMPM src/main.cpp src/utils.cpp src/Simulator.cpp)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(myMPM PRIVATE Eigen3::Eigen)

# 包含子目录，子目录里面也有 cmakeLists.txt
add_subdirectory(external/partio)
add_subdirectory(external/spdlog)

# 用在 add_executable() 之后
target_link_libraries(myMPM PRIVATE partio)
target_link_libraries(myMPM PRIVATE spdlog)
