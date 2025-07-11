set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(Boost REQUIRED)
message(STATUS "Boost include path : ${Boost_INCLUDE_DIR}")
message(STATUS "Boost library path : ${Boost_LIBRARY}")

set(CMAKE_CXX_FLAGS "${MCFF_WARNING_FLAGS} ${MCFF_OPTIMIZATION_FLAGS} ${MCFF_DEBUGGER_FLAGS}")
message(STATUS "CXX compiler flags are ${CMAKE_CXX_FLAGS}")

add_library(mcfflib STATIC
	ContactDetector.cpp
	Coordinate.cpp
	Input.cpp
	)
target_link_libraries(mcfflib PRIVATE
#	Boost::program_options)
	boost_program_options)
add_executable(mcff main.cpp)
set_target_properties(mcff PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
#target_include_directories(mcff PRIVATE Boost::program_options)
target_include_directories(mcff PRIVATE boost_program_options)
target_link_libraries(mcff PRIVATE
#	Boost::program_options
	boost_program_options
	mcfflib
	)

add_subdirectory(io)
add_subdirectory(util)

