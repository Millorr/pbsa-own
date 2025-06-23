#pragma once

#include "OpenGLRenderer.hpp"
#include "constants.hpp"

#include <OpenGLObjects.h>

#include <QPoint>
#include <QElapsedTimer>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>

class Prog2_3Simulation : public OpenGLRenderer
{
	Q_OBJECT

public:
	Prog2_3Simulation();
	~Prog2_3Simulation();

	void resize(int w, int h) override;
	void render() override;

	void mouseEvent(QMouseEvent * e) override;

	enum class BoundaryCondition
	{
		SINGLE_CORNER,
		TWO_CORNERS,
		DIAGONAL_CORNERS,
		THREE_CORNERS,
		ALL_CORNERS
	};

	Eigen::VectorXd step();

	void reset(
		double dt,
		BoundaryCondition boundaryCondition,
		double structKs, double structKd,
		double sheerKs, double sheerKd,
		double bendingKs, double bendingKd,
		int maxSimSteps
	);


private:
	double
		cameraAzimuth = constants::pi<double>/4,
		cameraElevation = constants::pi<double>/4;
	double cameraDistance = grid_size * dx; // Abstand der Kamera vom Ursprung

	bool rotateInteraction = false;
	bool firstStep = true;
	qint64 maxSimSteps = 10;

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
	double dt = 0.001;


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

	spring structureSpring = spring(100.0, 0.5, spring::Type::Structure); // Struktur-Federkonstante und Dämpfung
	spring sheerSpring = spring(50.0, 0.25, spring::Type::Sheer); // Scher-Federkonstante und Dämpfung
	spring bendingSpring = spring(10.0, 0.05, spring::Type::Bending); // Biege-Federkonstante und Dämpfung

	// Zustand
	Eigen::VectorXd positions;    // 3n Dimensionen (x,y,z für jeden Punkt)
	Eigen::VectorXd velocities;   // 3n Dimensionen
	Eigen::VectorXd forces;       // 3n Dimensionen

	// Hilfsfunktionen
	void computeForces();         // Berechnet alle Federkräfte
	void addSpringForce(int i1, int j1, int i2, int j2, const spring & spr);
	void midpointStep(); // Ein RK2-Schritt
	int getIndex(int i, int j) const;
	int getVectorIndex(int i, int j)const;
	int getVectorIndex(int particle_idx) const;
	void buildMesh();
	void buildCloth();
	void applyBoundaryConditions();
};
