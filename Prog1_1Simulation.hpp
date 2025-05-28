#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QPoint>
#include <QElapsedTimer>
#include <Eigen/Core>

class Prog1_1Simulation : public OpenGLRenderer
{
	Q_OBJECT

public:
	Prog1_1Simulation();
	~Prog1_1Simulation();

	Eigen::Vector3d  explicitEuler();
	Eigen::Vector3d  implicitEuler();
	Eigen::Vector3d  verlet();

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

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

private:
	double
		cameraAzimuth = constants::pi<double>,
		cameraElevation = constants::half_pi<double>;

	bool rotateInteraction = false;
	bool firstStep = true;

	QPointF lastPos;
	QElapsedTimer timer;
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

	gl::Texture
		earthTexture,
		moonTexture;

	GLsizei numIcosphereIndices = 0;

	// Simulation parameters
	Eigen::Vector3d r0 = Eigen::Vector3d(0.0,15.0,0.0), dr0 = Eigen::Vector3d(20.0, 0.0, 0.0);
	double g = 20.0, dt=0.01;
	Integration integration = Integration::Explicit;

	// Simulation state
	Eigen::Vector3d r = Eigen::Vector3d(0.0,15.0,0.0), dr = Eigen::Vector3d(20.0, 0.0, 0.0);
	Eigen::Vector3d r_prev;
};
