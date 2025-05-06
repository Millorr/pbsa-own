#pragma once

#include <Eigen/Core>

class Prog1_1Simulation
{
public:
	Prog1_1Simulation();
	~Prog1_1Simulation();

	enum class Integration
	{
		Explicit,
		Implicit,
		Verlet
	};

	Eigen::Vector3d step();
	void reset(
		Integration integration,
		Eigen::Vector3d r0,
		Eigen::Vector3d dr0,
		double g,
		double dt
	);
};
