#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QPoint>
#include <QElapsedTimer>
#include <QTimer>
#include <Eigen/Core>
#include <vector>

class Prog1_2Simulation : public OpenGLRenderer
{
	Q_OBJECT
public:
	Prog1_2Simulation();
	~Prog1_2Simulation();

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

	std::vector<Eigen::Vector3d> const & step();
	void reset(
		std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> const & bodies,
		double G,
		double dt
	);
private:
	double
		cameraAzimuth = constants::pi<double>,
		cameraElevation = constants::half_pi<double>;

	bool rotateInteraction = false;

	QPointF lastPos;
	QElapsedTimer timer;
	QTimer simTimer;
	quint64 lastTimeNS = 0;
	quint64 lastSimSteps = 0;

	int width = 0, height = 0;

	Eigen::Matrix4d
		projectionMatrix,
		viewMatrix;

	gl::Buffer
		icosphereVertexBuffer,
		icosphereIndexBuffer;

	gl::VertexArray
		icosphereVAO,
		skyboxVAO;

	gl::Program
		icosphereProgram,
		skyboxProgram;

	std::vector<gl::Texture> objectTextures;

	GLsizei numIcosphereIndices = 0;

	double G = 6.67430e-11, dt = 0.01;

	std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> bodies;
};
