#include "Prog2_2Simulation.hpp"

#include "helpers.hpp"

#include <QDebug>
#include <QMouseEvent>

#include <Eigen/Dense>
#include <Eigen/Geometry>

Prog2_2Simulation::Prog2_2Simulation()
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

Prog2_2Simulation::~Prog2_2Simulation() = default;


void Prog2_2Simulation::buildMesh()
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
			float x = dx * j - dx * (grid_size - 1) / 2;
			float y = dx * i - dx * (grid_size - 1) / 2;
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

void Prog2_2Simulation::buildCloth()
{

	// Vektoren initialisieren
	positions.resize(3 * grid_size * grid_size);
	velocities.resize(3 * grid_size * grid_size);

	velocities.setZero();  // Alle Geschwindigkeiten auf 0

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

void Prog2_2Simulation::resize(int w, int h)
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

void Prog2_2Simulation::render()
{
	auto currentTimeNS = this->timer.nsecsElapsed();
	auto deltaTimeNS = currentTimeNS - this->lastTimeNS;
	this->lastTimeNS = currentTimeNS;
	qint64 simSteps = (currentTimeNS / this->dt / 1e9) - this->lastSimSteps;
	this->lastSimSteps += simSteps;

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

		// step so oft ausführen, wie nötig, um die Simulation zu aktualisieren
		for(int i = 0; i < simSteps; i++)
		{
			this->step();
		}
		//qDebug() << "FPS: " << 1.0 / (deltaTimeNS / 1.0e09) << " simsteps: " << simSteps << " center: " << this->u(getIndex(grid_size/2, grid_size/2)) << " sum: " << this->u.sum();
		//qDebug() << "Kamera: " << cameraAzimuth << " eli: " << cameraElevation;

		//auto u_max = this->u.maxCoeff();

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

void Prog2_2Simulation::mouseEvent(QMouseEvent * e)
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

Eigen::VectorXd Prog2_2Simulation::step()
{
	if(firstStep)
	{
		lastTimeNS = timer.nsecsElapsed();
		firstStep = false;
		return positions;  // Ersten Frame überspringen
	}

	// do stuff
	applyBoundaryConditions();

	return positions;
}


void Prog2_2Simulation::applyBoundaryConditions()
{
	if(boundaryCondition == BoundaryCondition::SINGLE_CORNER)
	{
		positions.segment<3>(getVectorIndex(0, 0)) = Eigen::Vector3d(0.0, 0.0, 0.0);
		velocities.segment<3>(getVectorIndex(0, 0)).setZero();
	}
	else if(boundaryCondition == BoundaryCondition::BOTH_CORNERS)
	{
		positions.segment<3>(getVectorIndex(0, 0)) = Eigen::Vector3d(0.0, 0.0, 0.0);
		velocities.segment<3>(getVectorIndex(0, 0)).setZero();

		positions.segment<3>(getVectorIndex(32, 0)) = Eigen::Vector3d(32 * dx, 0.0, 0.0);
		velocities.segment<3>(getVectorIndex(32, 0)).setZero();
	}
}

int Prog2_2Simulation::getIndex(int i, int j) const
{
	// von 2d index auf 1d index
	return i * grid_size + j;
}

int Prog2_2Simulation::getVectorIndex(int i, int j) const
{
	return getIndex(i, j) * 3;  // Start-Index für x,y,z
}

int Prog2_2Simulation::getVectorIndex(int particle_idx) const
{
	return particle_idx * 3;  // Start-Index für x,y,z
}



void Prog2_2Simulation::reset(
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
	this->lastSimSteps = 0;
}
