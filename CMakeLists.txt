cmake_minimum_required (VERSION 3.15)

project (Mash VERSION 2.3)

set (CMAKE_STATIC_LIBRARY_PREFIX "")

string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE
         Release
         CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

# Compiler support.
include (CheckCXXCompilerFlag)

set (MASH_NATIVE_BUILD
     ON
     CACHE BOOL "Optimize build for current architecture.")
if ("${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Coverage")
    message (STATUS "${FontBold}Native build disabled due to Debug/Coverage build.${FontReset}")
elseif (MASH_NATIVE_BUILD)
    message (STATUS "${FontBold}Native build enabled. Built binaries will be optimized for this system.${FontReset}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else ()
    message (STATUS "${FontBold}Native build disabled. Detecting popcnt support.${FontReset}")
    check_cxx_compiler_flag ("-mpopcnt" MASH_HAS_POPCNT)
    if (MASH_HAS_POPCNT)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif ()
endif ()

check_cxx_compiler_flag ("-fopenmp" MASH_HAS_OPENMP)
if (MASH_HAS_OPENMP)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" MASH_HAS_OPENMP_SIMD)
if (MASH_HAS_OPENMP_SIMD)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
endif ()

check_cxx_compiler_flag ("-Wno-psabi" MASH_SUPPRESS_GCC4_ABI)
if (MASH_SUPPRESS_GCC4_ABI)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-psabi")
endif ()

check_cxx_compiler_flag ("-Wl,--wrap=memcpy" MASH_WRAP_MEMCPY)
if (MASH_WRAP_MEMCPY)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${CMAKE_CURRENT_LIST_DIR}/src/mash/memcpyLink.h -Wl,--wrap=memcpy")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -include ${CMAKE_CURRENT_LIST_DIR}/src/mash/memcpyLink.h -Wl,--wrap=memcpy")
endif ()


# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Require compression libs.
find_package (ZLIB QUIET REQUIRED)
find_package (BZip2 QUIET REQUIRED)

# Add SeqAn3 and Raptor.
set (RAPTOR_SUBMODULES_DIR "${CMAKE_CURRENT_LIST_DIR}/lib/raptor/lib")
set (SEQAN3_SUBMODULES_DIR "${RAPTOR_SUBMODULES_DIR}")
set (SEQAN3_CEREAL
     ON
     CACHE BOOL "Require cereal to be present.")
find_package (SeqAn3 QUIET REQUIRED HINTS ${RAPTOR_SUBMODULES_DIR}/seqan3/build_system)
add_subdirectory (lib/raptor/src EXCLUDE_FROM_ALL)

# Allow to include CMake scripts from seqan3.
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_CLONE_DIR}/test/cmake/")
# Allow to include CMake scripts from raptor.
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/lib/raptor/cmake/")

# # Define cmake configuration flags to configure and build external projects with the same flags as specified for
# # this project.
set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "--no-warn-unused-cli")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -w")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# Use ccache.
include (raptor_require_ccache)
raptor_require_ccache ()

# Dependency: xxHash
set (MASH_XXHASH_PATH "${PROJECT_BINARY_DIR}/lib/libxxhash.a")

if (EXISTS "${MASH_XXHASH_PATH}")
    add_custom_target (xxhash_project)
else ()
    set (MASH_XXHASH_ARGS "${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}")
    list (APPEND MASH_XXHASH_ARGS "-DXXHASH_BUILD_XXHSUM=OFF")
    list (APPEND MASH_XXHASH_ARGS "-DXXHASH_BUNDLED_MODE=OFF")
    list (APPEND MASH_XXHASH_ARGS "-DBUILD_SHARED_LIBS=OFF")
    list (APPEND MASH_XXHASH_ARGS "--no-warn-unused-cli")

    include (ExternalProject)
    ExternalProject_Add (
        xxhash_project
        PREFIX xxhash_project
        GIT_REPOSITORY https://github.com/Cyan4973/xxHash.git
        GIT_TAG "f9155bd4c57e2270a4ffbb176485e5d713de1c9b"
        CMAKE_ARGS "${MASH_XXHASH_ARGS}"
        SOURCE_SUBDIR cmake_unofficial
        BUILD_BYPRODUCTS "${MASH_XXHASH_PATH}"
        UPDATE_DISCONNECTED 1)
endif ()


file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")

add_library (xxhash STATIC IMPORTED)
add_dependencies (xxhash xxhash_project)
set_target_properties (xxhash PROPERTIES IMPORTED_LOCATION "${MASH_XXHASH_PATH}")
target_include_directories (xxhash INTERFACE "${PROJECT_BINARY_DIR}/include")

# Run Cap'n Proto on MinHash.
find_package (CapnProto)
include (cmake/CapnProtoMacros.cmake)
capnp_generate_cpp (src/mash/capnp/MinHash.capnp.c++ src/mash/capnp/MinHash.capnp.h src/mash/capnp/MinHash.capnp)

# Build mash.
add_library (libmash STATIC
             src/mash/capnp/MinHash.capnp.c++
             src/mash/Command.cpp
             src/mash/CommandBounds.cpp
             src/mash/CommandContain.cpp
             src/mash/CommandDistance.cpp
             src/mash/CommandFind.cpp
             src/mash/CommandIndex.cpp
             src/mash/CommandInfo.cpp
             src/mash/CommandList.cpp
             src/mash/CommandPaste.cpp
             src/mash/CommandScreen.cpp
             src/mash/CommandSketch.cpp
             src/mash/CommandTaxScreen.cpp
             src/mash/CommandTriangle.cpp
             src/mash/hash.cpp
             src/mash/HashList.cpp
             src/mash/HashPriorityQueue.cpp
             src/mash/HashSet.cpp
             src/mash/mash.cpp
             src/mash/MinHashHeap.cpp
             src/mash/MurmurHash3.cpp
             src/mash/Sketch.cpp
             src/mash/sketchParameterSetup.cpp
             src/mash/memcpyWrap.c)

target_link_libraries (libmash PUBLIC raptor_interface)
set_source_files_properties (src/mash/memcpyWrap.c PROPERTIES LANGUAGE C COMPILE_OPTIONS "-std=c17;-w")

target_link_libraries (libmash PUBLIC xxhash)

target_include_directories (libmash PUBLIC src)
target_include_directories (libmash PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_include_directories (libmash SYSTEM PUBLIC ${RAPTOR_SUBMODULES_DIR}/simde/simde)

target_compile_definitions (libmash PUBLIC "-DRAPTOR_HIBF_HAS_COUNT=1")

# Disable various warnings for libmash.
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-reorder")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-format")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-vla")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-char-subscripts")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-implicit-fallthrough")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-copy")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-catch-value")

# Add executable.
add_executable (mash src/mash/mash.cpp)
target_link_libraries (mash PUBLIC libmash)
