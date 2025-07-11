cmake_minimum_required(VERSION 3.10)
project(mcff)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD          17)

include_directories(${PROJECT_SOURCE_DIR})

include(CheckCXXCompilerFlag)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
find_package(Boost REQUIRED)
message(STATUS "Boost include path : ${Boost_INCLUDE_DIR}")
message(STATUS "Boost library path : ${Boost_LIBRARY}")

set(MCFF_WARNING_FLAGS "")
check_cxx_compiler_flag("-Wall" MCFF_COMPILER_HAS_WALL)
if(MCFF_COMPILER_HAS_WALL)
	set(MCFF_WARNING_FLAGS "${MCFF_WARNING_FLAGS} -Wall")
endif()

message(STATUS "warning flags are ${MCFF_WARNING_FLAGS}")

set(MCFF_OPTIMIZATION_FLAGS "")
option(DEBUG "compile in debug mode" OFF)
if(DEBUG)
	check_cxx_compiler_flag("-O0" MCFF_COMPILER_HAS_O0)
	if(MCFF_COMPILER_HAS_O0)
		set(MCFF_OPTIMIZATION_FLAGS "${MCFF_OPTIMIZATION_FLAGS} -O0")
	endif()
	check_cxx_compiler_flag("-g" MCFF_COMPILER_HAS_G)
	if(MCFF_COMPILER_HAS_G)
		set(MCFF_OPTIMIZATION_FLAGS "${MCFF_OPTIMIZATION_FLAGS} -g")
	endif()
else()
	check_cxx_compiler_flag("-Ofast" MCFF_COMPILER_HAS_OFAST)
	if(MCFF_COMPILER_HAS_OFAST)
		set(MCFF_OPTIMIZATION_FLAGS "${MCFF_OPTIMIZATION_FLAGS} -Ofast")
	endif()
endif()

message(STATUS "optimization flags are ${MCFF_OPTIMIZATION_FLAGS}")

set(MCFF_DEBUGGER_FLAGS "")
option(GDB "add flag for gdb" OFF)
if(GDB)
	check_cxx_compiler_flag("-ggdb" MCFF_COMPILER_HAS_GGDB)
	if(MCFF_COMPILER_HAS_GGDB)
		set(MCFF_DEBUGGER_FLAGS "${MCFF_DEBUGGER_FLAGS} -ggdb")
	endif()
endif()

message(STATUS "debugger flags are ${MCFF_DEBUGGER_FLAGS}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 8.0 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
		link_libraries(stdc++fs)
	endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
	link_libraries(stdc++fs)
endif()

add_subdirectory(src)
