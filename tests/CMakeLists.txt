project(SimulationFramework-Tests)

if(APPLE)
	set(TEST_MOC_HEADERS
		../OpenGLRenderer.hpp
	)

	qt6_wrap_cpp(TEST_MOC_SOURCES ${TEST_MOC_HEADERS})

	add_executable(
		${PROJECT_NAME}
		main.cpp
		../Prog1_1Simulation.cpp
		../Prog1_2Simulation.cpp
		../Prog1_3Simulation.cpp
		../helpers.cpp
		${TEST_MOC_SOURCES}
	)

	target_link_libraries(
		${PROJECT_NAME}
		PRIVATE
		Eigen3::Eigen
		Qt6::Core
		Qt6::Gui
		Qt6::Widgets
		Qt6::OpenGLWidgets
		OpenGLObjects
	)

	target_include_directories(
		${PROJECT_NAME}
		PRIVATE
		${CMAKE_SOURCE_DIR}
		${CMAKE_BINARY_DIR}
	)
else()
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
endif()
