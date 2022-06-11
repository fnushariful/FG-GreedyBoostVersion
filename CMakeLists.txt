#cmake_minimum_required(VERSION 3.21)
#project(FG-GreedyBoostVersion)

#set(CMAKE_CXX_STANDARD 17)

#set(CMAKE CXX FLAGS "03 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")
#add_executable(FG-GreedyBoostVersion main.cpp)



#find_package(CGAL)
#include_directories(/opt/homebrew/Cellar)
# Additional library search directories here
# Additional project files here main.cpp)
#target_link_libraries(FG-GreedyBoostVersion gmp)
#target_link_libraries(FG-GreedyBoostVersion mpfr)
#target_link_libraries(FG-GreedyBoostVersion boost_system)
#target_link_libraries(FG-GreedyBoostVersion boost_thread)
#target_link_libraries(FG-GreedyBoostVersion CGAL)


cmake_minimum_required(VERSION 3.21)
project(test)

include_directories(/opt/homebrew/Cellar)
include_directories(/opt/homebrew/Cellar)
include_directories(/opt/homebrew/Cellar/boost/1.76.0/include/boost)
#include_directories(/usr/local/Cellar/qt/6.2.2/include/)
#include_directories(/usr/local/Cellar/qt/6.2.2/QtWidgets/)
#include_directories(/usr/local/Cellar/qt/6.2.2/lib/cmake/)

link_directories(/opt/homebrew/lib/)
link_directories(/opt/homebrew/Cellar/)
#link_directories(/opt/homebrew/Cellar/boost/1.76.0/include/boost)
#link_directories(/usr/local/Cellar/qt/6.2.2/lib/cmake/)

find_package(CGAL)
find_package(Qt6Core REQUIRED)
find_package(Qt6Widgets REQUIRED)
find_package(Qt6Quick REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Xpreprocessor -fopenmp -std=c++17 -fexceptions -Wall -Wextra -O3 -stdlib=libc++ -Wno-unused-parameter -Wno-macro-redefined")
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/qt/6.2.2/lib/cmake/")
set(CMAKE_PREFIX_PATH "/opt/homebrew/Cellar/qt/6.2.2/QtWidgets/")

add_executable(test main.cpp GreedySpanner.h CGALComponents.h DelaunayTriangulation.h)

target_link_libraries(test libgmp.a)
target_link_libraries(test libmpfr.a)
target_link_libraries(test CGAL)
#target_link_libraries(FG-GreedyBoostVersion libomp.a)
target_link_libraries(test Qt6::Core)
target_link_libraries(test Qt6::Widgets)
target_link_libraries(test Qt6::Quick)
#target_link_libraries(FG-GreedyBoostVersion boost_system)
#target_link_libraries(FG-GreedyBoostVersion boost_thread)