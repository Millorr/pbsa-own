project(SimulationFramework-Tests)

add_executable(
	${PROJECT_NAME}
	main.cpp
	../Prog1_1Simulation.cpp
	../Prog1_2Simulation.cpp
	../Prog1_3Simulation.cpp
)

target_link_libraries(
	${PROJECT_NAME}
	PRIVATE
	Eigen3::Eigen
)
