#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QPoint>
#include <QElapsedTimer>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>

class Prog2_1Simulation : public OpenGLRenderer
{
	Q_OBJECT

public:
	Prog2_1Simulation();
	~Prog2_1Simulation();

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

	Eigen::VectorXd step();
	void reset(
		double dt,
		double a
	);
	void setTemperatureScale(
		double tempScale
	);


private:
	double
		cameraAzimuth = constants::pi<double>/4,
		cameraElevation = constants::pi<double>/4;

	bool rotateInteraction = false;
	bool firstStep = true;

	std::vector<float> vertices;
	std::vector<unsigned> indices;

	QPointF lastPos;
	QElapsedTimer timer;
	quint64 lastTimeNS = 0;
	quint64 lastSimSteps = 0;

	int width = 0, height = 0;

	Eigen::Matrix4d
		projectionMatrix,
		viewMatrix;

	gl::Buffer
		gridVertexBuffer,
		gridIndexBuffer;

	gl::VertexArray
		gridVAO,
		skyboxVAO;

	gl::Program
		gridProgram,
		skyboxProgram;

	gl::Texture
		gridTexture,
		skyboxTexture;
		

	GLsizei numGridIndices = 0;

	// Konstanten
	static constexpr int grid_size = 101;
	static constexpr double dx = 1.0;

	// Simulation Parameters
	double dt = 0.1;
	double a = 0.1;
	double tempScale = 5.0;

	// Zustand
	Eigen::VectorXd u;

	// Numerische Hilfsmittel
	Eigen::SparseMatrix<double> system_matrix;  // (I - a*dt*L)
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // Recommended for very sparse and not too large problems (e.g., 2D Poisson eq.) (ref: https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html)

	// Hilfsfunktionen
	int getIndex(int i, int j) const;
	void buildSystemMatrix();
	void setInitialConditions();
	Eigen::Vector3d calculateNormal(int i, int j) const;
};
