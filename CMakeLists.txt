cmake_minimum_required(VERSION 3.8)
project(f1tenth_planner)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)
find_package(yaml-cpp REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Python3_INCLUDE_DIRS}
)

# ===== matplotlib_cpp 라이브러리 =====
add_library(matplotlib_cpp INTERFACE)
target_include_directories(matplotlib_cpp INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(matplotlib_cpp INTERFACE cxx_std_11)
target_link_libraries(matplotlib_cpp INTERFACE
    Python3::Python
    Python3::Module
    yaml-cpp
)
if(Python3_NumPy_FOUND)
    target_link_libraries(matplotlib_cpp INTERFACE Python3::NumPy)
else()
    target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
endif()

# ===== main 실행 파일 =====
add_executable(main
    src/main.cpp
    src/visualization.cpp
    src/helper_func.cpp
    include/spline.cpp
)
target_link_libraries(main PRIVATE
    matplotlib_cpp
    Eigen3::Eigen
    yaml-cpp
)
set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ===== 설치 설정 =====
include(GNUInstallDirs)

set(PACKAGE_NAME matplotlib_cpp)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/${PACKAGE_NAME}/cmake)

install(TARGETS matplotlib_cpp EXPORT install_targets)
install(FILES
    "${PROJECT_SOURCE_DIR}/include/matplotlibcpp.h"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT install_targets
    FILE ${PACKAGE_NAME}Targets.cmake
    NAMESPACE ${PACKAGE_NAME}::
    DESTINATION ${INSTALL_CONFIGDIR}
)

# ===== 패키지 구성 파일 =====
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PACKAGE_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PACKAGE_NAME}Config.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}Config.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
