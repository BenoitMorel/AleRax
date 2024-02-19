
set(alerax_SOURCES ale.cpp 
  AleArguments.cpp
  AleEvaluator.cpp
  AleModelParameters.cpp
  AleOptimizer.cpp
  AleState.cpp
  Highways.cpp
  OptimizationClasses
  TrimFamilies.cpp
)
add_program(alerax "${alerax_SOURCES}")
target_link_libraries(alerax PRIVATE mpi-scheduler) 
install(TARGETS alerax DESTINATION bin)

