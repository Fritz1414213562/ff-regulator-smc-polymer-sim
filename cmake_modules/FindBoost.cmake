find_path(Boost_INCLUDE_DIR
	NAMES
		boost/program_options.hpp
	PATHS
		${Boost_ROOT}
	PATH_SUFFIXES
		include
	NO_DEFAULT_PATH
	)

find_path(Boost_LIBRARY
	NAMES
		libboost_program_options.a
		libboost_program_options.dylib
	PATHS
		${Boost_ROOT}
	PATH_SUFFIXES
		lib
	NO_DEFAULT_PATH
	)

#find_package(Boost REQUIRED COMPONENTS program_options)

mark_as_advanced(Boost_INCLUDE_DIR Boost_LIBRARY)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Boost
	REQUIRED_VARS
		Boost_LIBRARY
		Boost_INCLUDE_DIR
	)

if (Boost_FOUND AND NOT TARGET Boost::program_options)
	add_library(Boost::program_options UNKNOWN IMPORTED)
	set_target_properties(Boost::program_options PROPERTIES
		IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
		IMPORTED_LOCATION                 "${Boost_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES     "${Boost_INCLUDE_DIR}"
	)
endif()
