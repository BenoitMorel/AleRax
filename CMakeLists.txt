cmake_minimum_required(VERSION 3.0.1)

project(AleRaxProject
        VERSION 0.0)

set(ALERAX_BUILD_AS_LIBRARY OFF CACHE BOOL "Build AleRax as shared library (instead of stand-alone executable)")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

if(ALERAX_BUILD_AS_LIBRARY)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DALERAX_BUILD_AS_LIB=ON")
endif()

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++14 CXX_STANDARD_14_SUPPORT -g)
if(CXX_STANDARD_14_SUPPORT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
else()
    message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
      "Debug"
            CACHE STRING "Build type" FORCE
            )
          message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to Debug")
endif ()



add_subdirectory(ext)

if (NOT DISABLE_MPI)
  include(FindMPI)
  if(MPI_FOUND)
    set(WITH_MPI TRUE)
    add_definitions(-DWITH_MPI)
    message(STATUS "WITH MPI ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()


set(ALERAX_INCLUDE_DIRS ${GENERAXCORE_INCLUDE_DIRS}
   ${CMAKE_CURRENT_SOURCE_DIR}/ext
   ${MPI_INCLUDE_PATH}
   ${CORAX_INCLUDE_DIRS})
  

add_subdirectory(src)

