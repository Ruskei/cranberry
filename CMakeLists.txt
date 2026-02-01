cmake_minimum_required(VERSION 3.16)
project(cranberry VERSION 0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -march=native")
add_compile_options(-UNDEBUG)

find_package(VTK COMPONENTS
	CommonCore
	CommonDataModel
	FiltersCore
	IOXML
	ImagingMath
)

if (NOT VTK_FOUND)
	message(FATAL_ERROR "Unable to find VTK build folder.")
endif()

add_executable(cranberry
	src/main.cpp
	src/particles.cpp
	src/electrostatics.cpp
)

target_include_directories(cranberry
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cranberry
    PRIVATE
        ${VTK_LIBRARIES}
)

vtk_module_autoinit(
	TARGETS cranberry
	MODULES ${VTK_LIBRARIES}
)
