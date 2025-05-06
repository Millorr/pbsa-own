#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QElapsedTimer>
#include <QPoint>

#include <Eigen/Core>

class ExampleRenderer : public OpenGLRenderer
{
	Q_OBJECT

public:
	ExampleRenderer(QObject * parent);

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

private:
	double
		cameraAzimuth = constants::pi<double>,
		cameraElevation = constants::half_pi<double>;
	bool rotateInteraction = false;

	int width = 0, height = 0;

	QPointF lastPos;
	QElapsedTimer timer;
	quint64 lastTimeNS = 0;

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
		moonTexture,
		starsCubeMap;

	GLsizei numIcosphereIndices = 0;
};
