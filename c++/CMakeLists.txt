cmake_minimum_required(VERSION 3.5...3.19)
cmake_policy(SET CMP0048 NEW)
project(
  AI
  VERSION 0.0.1
  LANGUAGES "CXX" "C"
  DESCRIPTION "Artificial intelligence TSP solver")


find_package(Boost)


add_library(AI INTERFACE)
target_include_directories(AI INTERFACE core/)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS}) 
  target_link_libraries(AI ${Boost_LIBRARIES})
endif()


add_executable(main drivers/main.cpp)
target_link_libraries(main AI)

target_compile_features(main PRIVATE cxx_std_17)

