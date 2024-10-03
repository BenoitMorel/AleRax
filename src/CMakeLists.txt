
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

macro(add_program program_name program_sources)
  if (ALERAX_BUILD_AS_LIBRARY)
    add_library(${program_name} SHARED ${program_sources})
  else()
    add_executable(${program_name} ${program_sources})
  endif()
  target_sources(${program_name} PRIVATE ${program_sources})
  target_link_libraries(${program_name}
    PRIVATE
    generaxcore
    ${CORAX_LIB_TARGET}
    ${MPI_LIBRARIES}
    )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  if (ALERAX_BUILD_AS_LIBRARY AND NOT APPLE)
    target_link_libraries(${program_name} PRIVATE gcc_s gcc)
  endif()

  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name}
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name}
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )

  target_include_directories(${program_name}
    PUBLIC ${ALERAX_INCLUDE_DIRS}
    )
endmacro()

macro(add_program_corax program_name program_sources)
  if (ALERAX_BUILD_AS_LIBRARY)
    add_library(${program_name} SHARED ${program_sources})
  else()
    add_executable(${program_name} ${program_sources})
  endif()
  target_link_libraries(${program_name}
    generaxcore
    ${CORAX_LIB_TARGET}
    ${MPI_LIBRARIES}
    )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  if (ALERAX_BUILD_AS_LIBRARY AND NOT APPLE)
    target_link_libraries(${program_name} gcc_s gcc)
  endif()

  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name}
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name}
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )

  target_include_directories(${program_name}
    PUBLIC ${ALERAX_INCLUDE_DIRS}
    )
endmacro()



add_subdirectory(tools)
add_subdirectory(ale)
