
set(alerax_SOURCES ale.cpp
  AleArguments.cpp
  AleEvaluator.cpp
  AleOptimizer.cpp
  AleState.cpp
  Highways.cpp
  OptimizationClasses.cpp
  TrimFamilies.cpp
  )

add_program(alerax "${alerax_SOURCES}")
install(TARGETS alerax DESTINATION bin)

