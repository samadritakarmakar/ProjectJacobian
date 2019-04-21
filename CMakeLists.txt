cmake_minimum_required(VERSION 2.8)

project(ProjectJacobian)
find_package(Armadillo REQUIRED)
add_executable(${PROJECT_NAME} main.cpp f1.hpp f2.hpp Jacobian.hpp)
target_link_libraries(${PROJECT_NAME} armadillo)
