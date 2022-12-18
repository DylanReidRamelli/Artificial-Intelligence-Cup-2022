cmake_minimum_required(VERSION 3.5...3.19)
cmake_policy(SET CMP0048 NEW)
project(
  ACS
  VERSION 0.0.1
  LANGUAGES "C"
  DESCRIPTION "Artificial intelligence TSP solver")


add_library(ACS INTERFACE)
target_include_directories(ACS INTERFACE core/)
target_include_directories(ACS INTERFACE utils/)


add_executable(main drivers/main.c)
target_link_libraries(main ACS)
