#pragma once

#include <Eigen/Core>

#include <vector>

class Prog1_2Simulation
{
public:
	Prog1_2Simulation();
	~Prog1_2Simulation();

	std::vector<Eigen::Vector3d> const & step();
	void reset(
		std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> const & bodies,
		double G,
		double dt
	);
};
