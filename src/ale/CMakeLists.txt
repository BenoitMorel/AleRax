
set(alerax_SOURCES ale.cpp 
  AleArguments.cpp
  AleEvaluator.cpp
  AleOptimizer.cpp
  TrimFamilies.cpp
)
add_program(alerax "${alerax_SOURCES}")
target_link_libraries(alerax PRIVATE mpi-scheduler) 


