#include "Prog2_1Simulation.hpp"

#include "helpers.hpp"

#include <QDebug>
#include <QMouseEvent>

#include <Eigen/Dense>
#include <Eigen/Geometry>

Prog2_1Simulation::Prog2_1Simulation()
	: OpenGLRenderer{nullptr}
{


	{
		// create/bind a vertex array object to store bindings to array and element buffers
		glBindVertexArray(this->gridVAO.id());

		buildMesh();
		buildCloth();

		// create/bind an array buffer object and fill it with the vertices we computed (GL_DYNAMIC_DRAW, as it needs to be modified)
		glBindBuffer(GL_ARRAY_BUFFER, this->gridVertexBuffer.id());
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gridIndexBuffer.id());
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_DYNAMIC_DRAW);


		// Vertex Attribute für Position
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
		glEnableVertexAttribArray(0);

		// unbind the vertex array object as we are done modifying it
		glBindVertexArray(0);

		// remember the number of indices used (don't need the actual indices/vertices as these are on the GPU now)
		this->numGridIndices = static_cast<GLsizei>(indices.size());
	}

	{
		// create temporary shader objects of the correct types
		gl::Shader vertexShader{GL_VERTEX_SHADER};
		gl::Shader fragmentShader{GL_FRAGMENT_SHADER};

		std::vector<char> text;

		// load and compile vertex shader
		text = loadResource("shaders/clothShader.vert");
		vertexShader.compile(text.data(), static_cast<GLint>(text.size()));

		// load and compile fragment shader
		text = loadResource("shaders/clothShader.frag");
		fragmentShader.compile(text.data(), static_cast<GLint>(text.size()));

		// link shader program and check for errors. shader objects are no longer required
		if(!this->gridProgram.link(vertexShader, fragmentShader))
		{
			qDebug() << "Shader compilation failed:\n" << this->gridProgram.infoLog().get();
			std::abort();
		}

		// set which vertex attribute index is bound to position
		glBindAttribLocation(this->gridProgram.id(), 0, "position");
	}

	{
		// create temporary shader objects of the correct types
		gl::Shader vertexShader{GL_VERTEX_SHADER};
		gl::Shader fragmentShader{GL_FRAGMENT_SHADER};

		std::vector<char> text;

		// load and compile vertex shader
		text = loadResource("shaders/skybox.vert");
		vertexShader.compile(text.data(), static_cast<GLint>(text.size()));

		// load and compile fragment shader
		text = loadResource("shaders/skybox.frag");
		fragmentShader.compile(text.data(), static_cast<GLint>(text.size()));

		// link shader program and check for errors. shader objects are no longer required
		if(!this->skyboxProgram.link(vertexShader, fragmentShader))
		{
			qDebug() << "Shader compilation failed:\n" << this->skyboxProgram.infoLog().get();
			std::abort();
		}
	}
	{
		glBindTexture(GL_TEXTURE_2D, this->gridTexture.id());

		// 1x1 white pixel (RGBA)
		const uint8_t whitePixel[4] = {255, 255, 255, 255};

		// Upload the texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8_ALPHA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, whitePixel);

		// Set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		if(GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 1);
	}
	{
		glBindTexture(GL_TEXTURE_2D, this->skyboxTexture.id());

		// 1x1 white pixel (RGBA)
		const uint8_t whitePixel[4] = {0, 0, 0, 255};

		// Upload the texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8_ALPHA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, whitePixel);

		// Set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		if(GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 1);
	}

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.);

	glDisable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);


	this->timer.start();
}

Prog2_1Simulation::~Prog2_1Simulation() = default;


void Prog2_1Simulation::buildMesh()
{
	// sehr ineffizient, aber einfach zu verstehen
		// build the required vertex and index arrays
	vertices.clear();
	indices.clear();

	for(int i = 0; i < grid_size; ++i)
	{
		for(int j = 0; j < grid_size; ++j)
		{
			// Position
			float x = dx * j -dx * (grid_size - 1) / 2;
			float y = dx * i -dx * (grid_size - 1) / 2;
			float z = 0.0f;
			vertices.push_back(x);
			vertices.push_back(y);
			vertices.push_back(z);
		}
	}
	// Indices für zwei Dreiecke pro Zelle
	for(int i = 0; i < grid_size - 1; ++i)
	{
		for(int j = 0; j < grid_size - 1; ++j)
		{
			int idx = i * grid_size + j;
			indices.push_back(idx);
			indices.push_back(idx + 1);
			indices.push_back(idx + grid_size);

			indices.push_back(idx + 1);
			indices.push_back(idx + grid_size + 1);
			indices.push_back(idx + grid_size);
		}
	}
}

void Prog2_1Simulation::buildCloth()
{
	// Vektoren initialisieren
	positions.resize(3 * grid_size * grid_size);
	velocities.resize(3 * grid_size * grid_size);
	forces.resize(3 * grid_size * grid_size);

	velocities.setZero();
	forces.setZero();


	// Anfangspositionen setzen (undeformiert, parallel zum Boden)
	for(int i = 0; i < grid_size; ++i)
	{
		for(int j = 0; j < grid_size; ++j)
		{
			int idx = getVectorIndex(i, j);

			// Position in Metern
			positions[idx + 0] = j * dx;
			positions[idx + 1] = i * dx;
			positions[idx + 2] = 0.0f;
		}
	}
}

void Prog2_1Simulation::resize(int w, int h)
{
	this->width = w;
	this->height = h;
	// update projection matrix to account for (potentially) changed aspect ratio
	this->projectionMatrix = calculateInfinitePerspective(
		0.78539816339744831, // 45 degrees in radians
		static_cast<double>(w) / h,
		0.01 // near plane (chosen "at random")
	);
}

void Prog2_1Simulation::render()
{
	auto currentTimeNS = this->timer.nsecsElapsed();
	auto deltaTimeNS = currentTimeNS - this->lastTimeNS;
	this->lastTimeNS = currentTimeNS;
	int loc = -1;

	// clear color and depth buffer to black
	glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// recompute view matrix (camera position and direction) from azimuth/elevation
	{
		auto sa = std::sin(this->cameraAzimuth);
		auto ca = std::cos(this->cameraAzimuth);
		auto se = std::sin(this->cameraElevation);
		auto ce = std::cos(this->cameraElevation);
		auto distance = this->cameraDistance;
		this->viewMatrix = calculateLookAtMatrix(
			distance * Eigen::Vector3d{se * ca, se * sa, ce},
			{0, 0, 0},
			{0, 0, 1}
		);
	}

	// set the shader to be used
	glUseProgram(this->gridProgram.id());

	// set which texture unit is bound to colorTexture
	loc = glGetUniformLocation(this->gridProgram.id(), "colorTexture");
	glUniform1i(loc, 0);

	// set projection matrix
	loc = glGetUniformLocation(this->gridProgram.id(), "projection");
	glUniformMatrix4fv(loc, 1, GL_FALSE, this->projectionMatrix.cast<float>().eval().data());

	// bind the sphere VAO for drawing
	glBindVertexArray(this->gridVAO.id());

	{
		Eigen::Affine3d whiteTransform = Eigen::AlignedScaling3d{1.0, 1.0, 1.0} *Eigen::Affine3d::Identity();

		this->step();

		for(int i = 0; i < grid_size; ++i)
		{
			for(int j = 0; j < grid_size; ++j)
			{
				int idx = getVectorIndex(i, j);

				vertices[idx]		= positions[idx];
				vertices[idx + 1]	= positions[idx + 1];
				vertices[idx + 2]	= positions[idx + 2];
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, gridVertexBuffer.id());
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());


		// set modelView matrix
		loc = glGetUniformLocation(this->gridProgram.id(), "modelView");
		glUniformMatrix4fv(loc, 1, GL_FALSE, (this->viewMatrix * whiteTransform.matrix()).cast<float>().eval().data());

		// draw the sphere's triangles as indexed elements (using unsigned int indices), VAO already bound
		glDrawElements(GL_TRIANGLES, this->numGridIndices, GL_UNSIGNED_INT, nullptr);
	}


	// continuous animation, so we always want to redraw
	this->update();
}

void Prog2_1Simulation::mouseEvent(QMouseEvent * e)
{
	auto type = e->type();
	auto pos = e->position();

	// begin rotation interaction
	if(type == QEvent::MouseButtonPress && e->button() == Qt::LeftButton)
	{
		this->lastPos = pos;
		this->rotateInteraction = true;
		return;
	}

	// end rotation interaction
	if(type == QEvent::MouseButtonRelease && e->button() == Qt::LeftButton)
	{
		this->rotateInteraction = false;
		return;
	}

	// perform rotation interaction
	if(this->rotateInteraction)
	{
		auto delta = pos - this->lastPos;
		this->lastPos = pos;

		// scale rotation depending on window diagonal
		auto scale = constants::two_pi<double> / Eigen::Vector2d{this->width, this->height}.norm();

		// modify azimuth and elevation by change in mouse position
		cameraAzimuth -= scale * delta.x();
		cameraAzimuth = std::fmod(cameraAzimuth, constants::two_pi<double>);
		cameraElevation -= scale * delta.y();

		// limit elevation so up is never collinear with camera view direction
		cameraElevation = std::fmax(std::fmin(cameraElevation, constants::pi<double> -0.01), 0.01);

		// tell widget to update itself to account for changed position
		this->update();
	}
}

Eigen::VectorXd Prog2_1Simulation::step()
{
	if(firstStep)
	{
		lastTimeNS = timer.nsecsElapsed();
		firstStep = false;
		return positions;  // Ersten Frame überspringen
	}

	midpointStep();
	applyBoundaryConditions();

	return positions;
}

void Prog2_1Simulation::computeForces()
{
	forces.setZero();  // Alle Kräfte zurücksetzen

	// Gravitationskraft für alle Punkte
	for(int i = 0; i < grid_size * grid_size; ++i)
	{
		forces[3 * i + 2] -= m * g;
	}

	// Federkräfte berechnen
	for(int i = 0; i < grid_size; ++i)
	{
		for(int j = 0; j < grid_size; ++j)
		{
			int idx = getIndex(i, j);

			// Strukturfedern (horizontal und vertikal)
			addSpringForce(i, j, i + 1, j, structureSpring);     // oben
			addSpringForce(i, j, i, j + 1, structureSpring);     // rechts

			// Scherfedern (diagonal)
			addSpringForce(i, j, i + 1, j + 1, sheerSpring);       // Diagonal rechts-oben
			addSpringForce(i, j, i + 1, j - 1, sheerSpring);       // Diagonal links-oben

			// Biegefedern (2er-Abstand)
			addSpringForce(i, j, i + 2, j, bendingSpring);       // 2 nach oben
			addSpringForce(i, j, i, j + 2, bendingSpring);       // 2 nach rechts
		}
	}
}

void Prog2_1Simulation::addSpringForce(int i1, int j1, int i2, int j2, const spring & spr)
{
	// Grenzen prüfen
	if(i2 < 0 || i2 >= grid_size || j2 < 0 || j2 >= grid_size)
	{
		return;
	}

	int idx1 = getVectorIndex(i1, j1);
	int idx2 = getVectorIndex(i2, j2);

	// Positionen und Geschwindigkeiten extrahieren
	Eigen::Vector3d pos1 = positions.segment<3>(idx1);
	Eigen::Vector3d pos2 = positions.segment<3>(idx2);
	Eigen::Vector3d vel1 = velocities.segment<3>(idx1);
	Eigen::Vector3d vel2 = velocities.segment<3>(idx2);

	// Relative Vektoren
	Eigen::Vector3d x_ij = pos2 - pos1;
	Eigen::Vector3d v_ij = vel2 - vel1;

	double x_ij_norm = x_ij.norm();

	if(x_ij_norm < 1e-12) return;  // Schutz vor Division durch 0

	Eigen::Vector3d x_ij_hat = x_ij / x_ij_norm;

	// Ruhelänge der Feder
	double rest_length = dx;
	if(spr.type == spring::Type::Structure)
	{
		rest_length = dx;  // Strukturfedern haben eine Ruhelänge von dx
	}
	else if(spr.type == spring::Type::Sheer)
	{
		rest_length = dx * sqrt(2.0); // Scherfedern haben eine Ruhelänge von dx * sqrt(2)
	}
	else if(spr.type == spring::Type::Bending)
	{
		rest_length = 2.0 * dx; // Biegefedern haben eine Ruhelänge von 2 * dx
	}

	// Federkraft berechnen (Gleichung 2 aus der Aufgabe)
	double spring_force_magnitude = spr.ks * (x_ij_norm - rest_length);
	double damping_force_magnitude = spr.kd * v_ij.dot(x_ij_hat);

	double total_force_magnitude = spring_force_magnitude + damping_force_magnitude;

	Eigen::Vector3d force = total_force_magnitude * x_ij_hat;

	// Kräfte auf beide Punkte anwenden (Newton's 3. Gesetz)
	forces.segment<3>(idx1) += force;
	forces.segment<3>(idx2) -= force;
}

void Prog2_1Simulation::midpointStep()
{
	// k1 = f(t, q) berechnen
	computeForces();

	Eigen::VectorXd k1_vel = forces / m;  // a = F/m
	Eigen::VectorXd k1_pos = velocities;

	// Zwischenschritt: q + dt/2 * k1
	Eigen::VectorXd temp_positions = positions + 0.5 * dt * k1_pos;
	Eigen::VectorXd temp_velocities = velocities + 0.5 * dt * k1_vel;

	// Temporären Zustand setzen
	Eigen::VectorXd original_positions = positions;
	Eigen::VectorXd original_velocities = velocities;

	positions = temp_positions;
	velocities = temp_velocities;

	// k2 = f(t + dt/2, q + dt/2 * k1) berechnen
	computeForces();

	Eigen::VectorXd k2_vel = forces / m;
	Eigen::VectorXd k2_pos = temp_velocities;

	// Originalen Zustand wiederherstellen
	positions = original_positions;
	velocities = original_velocities;

	// Finales Update mit k2
	velocities += dt * k2_vel;
	positions += dt * k2_pos;
}

void Prog2_1Simulation::applyBoundaryConditions()
{
	// Anfangs Ecke (0,0) immer fixieren
	positions.segment<3>(getVectorIndex(0, 0)) = Eigen::Vector3d(0.0, 0.0, 0.0);
	velocities.segment<3>(getVectorIndex(0, 0)).setZero();

	// Ecke nebenan (grid_size - 1, 0) fixieren
	if(boundaryCondition == BoundaryCondition::TWO_CORNERS || boundaryCondition == BoundaryCondition::THREE_CORNERS || boundaryCondition == BoundaryCondition::ALL_CORNERS)
	{
		positions.segment<3>(getVectorIndex(grid_size - 1, 0)) = Eigen::Vector3d((grid_size - 1) * dx, 0.0, 0.0);
		velocities.segment<3>(getVectorIndex(grid_size - 1, 0)).setZero();
	}

	// Diagonale Ecke (grid_size - 1, grid_size - 1) fixieren
	if(boundaryCondition == BoundaryCondition::DIAGONAL_CORNERS || boundaryCondition == BoundaryCondition::ALL_CORNERS)
	{
		positions.segment<3>(getVectorIndex(grid_size - 1, grid_size - 1)) = Eigen::Vector3d((grid_size - 1) * dx, (grid_size - 1) * dx, 0.0);
		velocities.segment<3>(getVectorIndex(grid_size - 1, grid_size - 1)).setZero();
	}

	// Letze Ecke (0, grid_size - 1) fixieren
	if(boundaryCondition == BoundaryCondition::ALL_CORNERS || boundaryCondition == BoundaryCondition::THREE_CORNERS)
	{
		positions.segment<3>(getVectorIndex(0, grid_size - 1)) = Eigen::Vector3d(0.0, (grid_size - 1) * dx, 0.0);
		velocities.segment<3>(getVectorIndex(0, grid_size - 1)).setZero();
	}
}

int Prog2_1Simulation::getIndex(int i, int j) const
{
	// von 2d index auf 1d index
	return i * grid_size + j;
}

int Prog2_1Simulation::getVectorIndex(int i, int j) const
{
	return getIndex(i, j) * 3;  // Start-Index für x,y,z
}

int Prog2_1Simulation::getVectorIndex(int particle_idx) const
{
	return particle_idx * 3;  // Start-Index für x,y,z
}

void Prog2_1Simulation::reset(
	double dt0Param,
	BoundaryCondition boundaryConditionParam,
	double structKsParam, double structKdParam,
	double sheerKsParam, double sheerKdParam,
	double bendingKsParam, double bendingKdParam
)
{
	this->dt = dt0Param;
	this->boundaryCondition = boundaryConditionParam;
	this->structureSpring = spring(structKsParam, structKdParam);
	this->sheerSpring = spring(sheerKsParam, sheerKdParam);
	this->bendingSpring = spring(bendingKsParam, bendingKdParam);

	qDebug() << "Reset simulation with parameters:"
		<< "dt0:" << dt
		<< "boundaryCondition:" << static_cast<int>(boundaryConditionParam)
		<< "structKs:" << structKsParam
		<< "structKd:" << structKdParam
		<< "sheerKs:" << sheerKsParam
		<< "sheerKd:" << sheerKdParam
		<< "bendingKs:" << bendingKsParam
		<< "bendingKd:" << bendingKdParam;


	buildCloth();

	this->lastTimeNS = this->timer.nsecsElapsed();
	this->timer.restart();
}
