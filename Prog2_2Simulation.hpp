#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QPoint>
#include <QElapsedTimer>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>

class Prog2_2Simulation : public OpenGLRenderer
{
	Q_OBJECT

public:
	Prog2_2Simulation();
	~Prog2_2Simulation();

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

	enum class BoundaryCondition
	{
		SINGLE_CORNER,    // Nur (0,0) fixiert
		BOTH_CORNERS      // (0,0) und (32,0) fixiert
	};

	Eigen::VectorXd step();

	void reset(
		double dt,
		BoundaryCondition boundaryCondition,
		double structKs, double structKd,
		double sheerKs, double sheerKd,
		double bendingKs, double bendingKd
	);


private:
	double
		cameraAzimuth = constants::pi<double> / 4,
		cameraElevation = constants::pi<double> / 4;
	double cameraDistance = 32.0; // Abstand der Kamera vom Ursprung

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
	static constexpr int grid_size = 33;
	static constexpr double dx = 0.10; // 10cm in m
	static constexpr double g = 10.0;
	static constexpr double m = 0.03;  // 30g in kg

	// Simulation Parameters
	double dt = 0.1;


	struct spring
	{
		enum class Type
		{
			Structure, Sheer, Bending
		};

		double ks; // Federkonstante
		double kd; // Dämpfungskonstante
		Type type;

		spring(double ks_ = 0.0, double kd_ = 0.0, Type t = Type::Structure)
			: ks(ks_), kd(kd_), type(t)
		{}
	};

	BoundaryCondition boundaryCondition = BoundaryCondition::SINGLE_CORNER;

	spring structureSpring = spring(2000.0, 5.0, spring::Type::Structure); // Struktur-Federkonstante und Dämpfung
	spring sheerSpring = spring(1000.0, 2.5, spring::Type::Sheer); // Scher-Federkonstante und Dämpfung
	spring bendingSpring = spring(200.0, 1.0, spring::Type::Bending); // Biege-Federkonstante und Dämpfung

	// Zustand
	Eigen::VectorXd positions;    // 3*N Einträge
	Eigen::VectorXd velocities;   // 3*N Einträge


	// Numerische Hilfsmittel
	Eigen::SparseMatrix<double> system_matrix;  // (I - a*dt*L)
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver; // Recommended for very sparse and not too large problems (e.g., 2D Poisson eq.) (ref: https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html)

	// Hilfsfunktionen
	int getIndex(int i, int j) const;
	int getVectorIndex(int i, int j)const;
	int getVectorIndex(int particle_idx) const;

	void buildSystemMatrix();
	void setInitialConditions();

	void buildMesh();
	void buildCloth();
	void applyBoundaryConditions();
	void setInitialPositions();
};
