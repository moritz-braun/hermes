add_subdirectory(quadrature)
add_subdirectory(bubbles)
add_subdirectory(mesh)
add_subdirectory(tutorial)
add_subdirectory(benchmarks)
add_subdirectory(examples)
add_subdirectory(adaptivity)
add_subdirectory(view)
add_subdirectory(shapeset)
add_subdirectory(integrals)

# Additional definitions for tests.
add_definitions(-DH2D_REPORT_ALL -DH2D_TEST)
