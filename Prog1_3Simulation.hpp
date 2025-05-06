#pragma once

#include <Eigen/Core>

class Prog1_3Simulation
{
public:
	Prog1_3Simulation();
	~Prog1_3Simulation();

	Eigen::VectorXd step();
	void reset(
		double dt,
		double a
	);
};
