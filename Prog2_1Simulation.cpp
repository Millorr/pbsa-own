#include "Prog2_1Simulation.hpp"

#include "helpers.hpp"

#include <QDebug>
#include <QMouseEvent>

#include <Eigen/Dense>
#include <Eigen/Geometry>

Prog2_1Simulation::Prog2_1Simulation()
	: OpenGLRenderer{nullptr}
{

	buildSystemMatrix();
	this->solver.compute(this->system_matrix);
	setInitialConditions();

	{
		// create/bind a vertex array object to store bindings to array and element buffers
		glBindVertexArray(this->gridVAO.id());

		// sehr ineffizient, aber einfach zu verstehen
		// build the required vertex and index arrays
		vertices.clear();
		indices.clear();
		for(int i = 0; i < grid_size; ++i)
		{
			for(int j = 0; j < grid_size; ++j)
			{
				// Position
				float x = -1.0f + 2.0f * j / (grid_size - 1);
				float y = -1.0f + 2.0f * i / (grid_size - 1);
				float z = u(getIndex(i,j));
				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);

				// Normals
				Eigen::Vector3d normal = calculateNormal(i, j);

				vertices.push_back(normal.x());
				vertices.push_back(normal.y());  
				vertices.push_back(normal.z());

				// Farben (blau bis rot)
				float value = std::clamp(static_cast<float>(u(getIndex(i, j))), 0.0f, 1.0f);
				float r = value;
				float g = 0.0f;
				float b = 1.0f - value;
				vertices.push_back(r);
				vertices.push_back(g);
				vertices.push_back(b);
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

		// create/bind an array buffer object and fill it with the vertices we computed (GL_DYNAMIC_DRAW, as it needs to be modified)
		glBindBuffer(GL_ARRAY_BUFFER, this->gridVertexBuffer.id());
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_DYNAMIC_DRAW);


		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->gridIndexBuffer.id());
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);


		// Vertex Attribute für Position
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)0);
		glEnableVertexAttribArray(0);

		// Vertex Attribute für Normalen
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)(3 * sizeof(float)));
		glEnableVertexAttribArray(1);

		// Vertex Attribute für Farbe
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)(6 * sizeof(float)));
		glEnableVertexAttribArray(2);

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
		text = loadResource("shaders/lambertShader.vert");
		vertexShader.compile(text.data(), static_cast<GLint>(text.size()));

		// load and compile fragment shader
		text = loadResource("shaders/lambertShader.frag");
		fragmentShader.compile(text.data(), static_cast<GLint>(text.size()));

		// link shader program and check for errors. shader objects are no longer required
		if(!this->gridProgram.link(vertexShader, fragmentShader))
		{
			qDebug() << "Shader compilation failed:\n" << this->gridProgram.infoLog().get();
			std::abort();
		}

		// set which vertex attribute index is bound to position
		glBindAttribLocation(this->gridProgram.id(), 0, "position");
		glBindAttribLocation(this->gridProgram.id(), 1, "normal");
		glBindAttribLocation(this->gridProgram.id(), 2, "color");
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

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);


	this->timer.start();
}

Prog2_1Simulation::~Prog2_1Simulation() = default;


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
		auto distance = 4.;
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
		Eigen::Affine3d whiteTransform = Eigen::AlignedScaling3d{1.0, 1.0, this->tempScale} *Eigen::Affine3d::Identity();

		// step so oft ausführen, wie nötig, um die Simulation zu aktualisieren
		for(int i = 0; i < simSteps; i++)
		{
			this->step();
		}
		//qDebug() << "FPS: " << 1.0 / (deltaTimeNS / 1.0e09) << " simsteps: " << simSteps << " center: " << this->u(getIndex(grid_size/2, grid_size/2)) << " sum: " << this->u.sum();
		//qDebug() << "Kamera: " << cameraAzimuth << " eli: " << cameraElevation;

		auto u_max = this->u.maxCoeff();

		for(int i = 0; i < grid_size; ++i)
		{
			for(int j = 0; j < grid_size; ++j)
			{
				int idx = (i * grid_size + j) * 9 + 2; // Z-Komponente
				// Aktualisiere die Z-Komponente der vertecies
				vertices[idx] = u(getIndex(i, j));

				// Aktualisiere die Normalen
				Eigen::Vector3d normal = calculateNormal(i, j);

				vertices[idx + 1] = normal.x();
				vertices[idx + 2] = normal.y();
				vertices[idx + 3] = normal.z();

				// Farbe neu berechnen (blau bis rot)
				float value = std::clamp(static_cast<float>(u(getIndex(i, j)) / u_max), 0.0f, 1.0f);
				float r = value;
				float g = 0.0f;
				float b = 1.0f - value;
				vertices[idx + 4] = r;
				vertices[idx + 5] = g;
				vertices[idx + 6] = b;
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
	this->u = solver.solve(this->u);
	return this->u;
}

int Prog2_1Simulation::getIndex(int i, int j) const
{
	// von 2d index auf 1d index
	return i * grid_size + j;
}

void Prog2_1Simulation::buildSystemMatrix()
{
	const int N = grid_size * grid_size;

	// Beginne mit Identitätsmatrix I
	// N * N Matrix, da wir GRID_SIZE x GRID_SIZE Knoten haben
	// Und wir abbilden wollen wie viel Einfluss jeder Knoten auf jeden anderen hat
	system_matrix = Eigen::SparseMatrix<double>(N, N);
	system_matrix.setIdentity();

	// Faktor für finite Differenzen: a*dt/dx²
	const double factor = a * dt; // default: 0.0000001

	// Durchlaufe alle Gitterpunkte
	for(int i = 0; i < grid_size; i++)
	{
		for(int j = 0; j < grid_size; j++)
		{
			int idx = getIndex(i, j);

			// Randbedingungen: u = 0 am gesamten Rand
			if(i == 0 || i == grid_size - 1 || j == 0 || j == grid_size - 1)
			{
				// Für Randknoten: Matrix-Zeile bleibt [0...0 1 0...0]
				// Das erzwingt u[idx] = 0 (Dirichlet-Randbedingung)
				continue;
			}


			// (I - dt*a*L) * u_t+1 = u_t
			// Hauptdiagonale: 1 + 4*factor (vom -4 im Laplace-Operator)
            system_matrix.coeffRef(idx, idx) = 1.0 + 4.0 * factor;

			// 5-Punkt-Stencil: Nachbarn mit -factor
			// Refernce zu den Nachbarn setzen (-1, 0), (1, 0), (0, -1), (0, 1)
			// 9 Punkte geht auch, aber weiß nicht so ganz, ob das richtig wäre, da es dann über dx von 1 hinaus gehen würde.
			system_matrix.coeffRef(idx, getIndex(i - 1, j)) = -factor;  // Links
			system_matrix.coeffRef(idx, getIndex(i + 1, j)) = -factor;  // Rechts
			system_matrix.coeffRef(idx, getIndex(i, j - 1)) = -factor;  // Unten  
			system_matrix.coeffRef(idx, getIndex(i, j + 1)) = -factor;  // Oben
		}
	}
}

void Prog2_1Simulation::setInitialConditions()
{
	this->u = Eigen::VectorXd(grid_size*grid_size);
	this->u.setZero(); // Setze alle Knoten auf 0.0 (Anfangstemperatur)

	int mid = grid_size / 2;
	int centerIndex = getIndex(mid, mid);
	// Setze den zentralen Knoten auf 1.0 (Anfangstemperatur)
	this->u(centerIndex) = 1.0;
}

Eigen::Vector3d Prog2_1Simulation::calculateNormal(int i, int j) const
{
	// Randbehandlung: Normale nach oben
	if(i <= 0 || i >= grid_size - 1 || j <= 0 || j >= grid_size - 1)
		return Eigen::Vector3d(0, 0, 1);

	// Zentrale Differenzen
	float dzdx = static_cast<float>(u(getIndex(i, j + 1)) - u(getIndex(i, j - 1))) / (2.0f * dx);
	float dzdy = static_cast<float>(u(getIndex(i + 1, j)) - u(getIndex(i - 1, j))) / (2.0f * dx);

	// Normale: (-dz/dx, -dz/dy, 1)
	Eigen::Vector3d n = {-dzdx, -dzdy, 1.0f};
	n.normalize();
	return n;
}

void Prog2_1Simulation::reset(
	double dt0Param,
	double a0Param
)
{
	this->dt = dt0Param;
	this->a = a0Param;


	buildSystemMatrix();
	this->solver.compute(this->system_matrix);
	setInitialConditions();

	qDebug() << "Reset simulation with parameters:"
		<< "dt0:" << dt
		<< "a0:" << a;

	this->lastTimeNS = this->timer.nsecsElapsed();
	this->timer.restart();
	this->lastSimSteps = 0;
}

void Prog2_1Simulation::setTemperatureScale(double tempScaleParam)
{
	this->tempScale = tempScaleParam;
	// qDebug() << "Set temperature scale to" << tempScaleParam;
}
