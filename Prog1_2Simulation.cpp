#include "Prog1_2Simulation.hpp"

#include "helpers.hpp"

#include <QImage>
#include <QDebug>
#include <QMouseEvent>
#include <QTimer>

#include <Eigen/Dense>
#include <Eigen/Geometry>



Prog1_2Simulation::Prog1_2Simulation()
	: OpenGLRenderer{nullptr}
{
	{
		// create/bind a vertex array object to store bindings to array and element buffers
		glBindVertexArray(this->icosphereVAO.id());

		// build the required vertex and index arrays
		std::vector<float> vertices(std::begin(icosahedronVertices), std::end(icosahedronVertices));
		std::vector<unsigned> indices(std::begin(icosahedronIndices), std::end(icosahedronIndices));
		for(auto k = 4; k--;) {
			subdivideIcosphere(vertices, indices);
		}

		// create/bind an array buffer object and fill it with the vertices we computed (GL_STATIC_DRAW, as this remains unmodified)
		glBindBuffer(GL_ARRAY_BUFFER, this->icosphereVertexBuffer.id());
		glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(vertices[0]) * vertices.size()), vertices.data(), GL_STATIC_DRAW);

		// set vertex attribute 0 to use the current array buffer
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

		// create/bind an element (index) buffer and fill it with the indices we computed (GL_STATIC_DRAW, as this remains unmodified)
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, this->icosphereIndexBuffer.id());
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(sizeof(indices[0]) * indices.size()), indices.data(), GL_STATIC_DRAW);

		// unbind the vertex array object as we are done modifying it
		glBindVertexArray(0);

		// remember the number of indices used (don't need the actual indices/vertices as these are on the GPU now)
		this->numIcosphereIndices = static_cast<GLsizei>(indices.size());
	}

	{
		// create temporary shader objects of the correct types
		gl::Shader vertexShader{GL_VERTEX_SHADER};
		gl::Shader fragmentShader{GL_FRAGMENT_SHADER};

		std::vector<char> text;

		// load and compile vertex shader
		text = loadResource("shaders/icosphere.vert");
		vertexShader.compile(text.data(), static_cast<GLint>(text.size()));

		// load and compile fragment shader
		text = loadResource("shaders/icosphere.frag");
		fragmentShader.compile(text.data(), static_cast<GLint>(text.size()));

		// link shader program and check for errors. shader objects are no longer required
		if(!this->icosphereProgram.link(vertexShader, fragmentShader))
		{
			qDebug() << "Shader compilation failed:\n" << this->icosphereProgram.infoLog().get();
			std::abort();
		}

		// set which vertex attribute index is bound to position
		glBindAttribLocation(this->icosphereProgram.id(), 0, "position");
	}

	{
		// Set bodies to default values
		//this->bodies = {
		//	{Eigen::Vector3d{15.0, 3.0, -3.0}, Eigen::Vector3d{2.0, -1.0, 0.0}, 1.0},
		//	{Eigen::Vector3d{11.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 0.0, 0.5}, 1.0},
		//	{Eigen::Vector3d{17.0, 6.0, 4.0}, Eigen::Vector3d{1.0, 0.0, 0.0}, 1.0}
		//};
		//this->reset(this->bodies, 6.67430e-11, 0.01);

		//values for a stable 3 body system ->Lagrange Dreieck
		std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> bodies_ = {
			// r0, v0, m
			{ Eigen::Vector3d{0.0, 1.0, 0.0},         Eigen::Vector3d{-1.224745, 0.0, 0.0},         1.0 },
			{ Eigen::Vector3d{0.8660254, -0.5, 0.0},  Eigen::Vector3d{0.612372, 1.06066, 0.0},      1.0 },
			{ Eigen::Vector3d{-0.8660254, -0.5, 0.0}, Eigen::Vector3d{0.612372, -1.06066, 0.0},     1.0 }
		};
		double G_ = 1.0;
		double dt_ = 0.01;
		this->reset(bodies_, G_, dt_);
		//this->reset(this->bodies, 6.67430e-11, 0.01);
	}

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	this->timer.start();
}

Prog1_2Simulation::~Prog1_2Simulation() = default;

void Prog1_2Simulation::resize(int w, int h)
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

void Prog1_2Simulation::render(){
	auto currentTimeNS = this->timer.nsecsElapsed();
	auto deltaTimeNS = currentTimeNS - this->lastTimeNS;
	this->lastTimeNS = currentTimeNS;
	qint64 simSteps = (currentTimeNS/this->dt/ 1e9) - this->lastSimSteps;
	this->lastSimSteps += simSteps;

	int loc = -1;

	std::vector<Eigen::Vector3d> positions;
	for(int i = 0; i < simSteps; i++)
	{
		positions = this->step();
	}
	if (positions.empty())
	{
		for (const auto & body : this->bodies)
		{
			Eigen::Vector3d position = std::get<0>(body);
			positions.push_back(position);
		}
	}

	// clear color and depth buffer to black
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
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
	glUseProgram(this->icosphereProgram.id());

	// set which texture unit is bound to colorTexture
	loc = glGetUniformLocation(this->icosphereProgram.id(), "colorTexture");
	glUniform1i(loc, 0);

	// set projection matrix
	loc = glGetUniformLocation(this->icosphereProgram.id(), "projection");
	glUniformMatrix4fv(loc, 1, GL_FALSE, this->projectionMatrix.cast<float>().eval().data());

	// bind the sphere VAO for drawing
	glBindVertexArray(this->icosphereVAO.id());

	qDebug() << "Simulation step took" << deltaTimeNS / 1e6 << "ms, simSteps:" << simSteps;
	for (unsigned long i = 0; i < this->bodies.size(); ++i)
	{
		Eigen::Vector3d position = positions[i];
		Eigen::Affine3d bodyTransform = Eigen::Translation3d{position.x(), position.y(), position.z()} * Eigen::Affine3d::Identity();

		// bind texture to unit 0
		glActiveTexture(GL_TEXTURE0 + 0);
		glBindTexture(GL_TEXTURE_2D, this->objectTextures[i].id());

		// set modelView matrix
		loc = glGetUniformLocation(this->icosphereProgram.id(), "modelView");
		glUniformMatrix4fv(loc, 1, GL_FALSE, (this->viewMatrix * bodyTransform.matrix()).cast<float>().eval().data());

		// draw the sphere's triangles as indexed elements (using unsigned int indices), VAO already bound
		glDrawElements(GL_TRIANGLES, this->numIcosphereIndices, GL_UNSIGNED_INT, nullptr);
	}

	this->update();
}

std::vector<Eigen::Vector3d> const & Prog1_2Simulation::step()
{
	// Placeholder
	//static std::vector<Eigen::Vector3d> positions;
	//positions.clear();
	//for (const auto & body : this->bodies)
	//{
	//	Eigen::Vector3d position = std::get<0>(body);
	//	positions.push_back(position);
	//}

	static std::vector<Eigen::Vector3d> positions;
	positions.clear();

	// Kopiere aktuelle Zustände
	std::vector<Eigen::Vector3d> newPositions;
	std::vector<Eigen::Vector3d> newVelocities;

	for(const auto & body : this->bodies)
	{
		newPositions.push_back(std::get<0>(body));
		newVelocities.push_back(std::get<1>(body));
	}

	// Berechne Kräfte und aktualisiere Positionen & Geschwindigkeiten
	for(size_t i = 0; i < this->bodies.size(); ++i)
	{
		Eigen::Vector3d force = Eigen::Vector3d::Zero();
		Eigen::Vector3d xi = std::get<0>(this->bodies[i]);
		double mi = std::get<2>(this->bodies[i]);

		for(size_t j = 0; j < this->bodies.size(); ++j)
		{
			if(i == j) continue;
			Eigen::Vector3d xj = std::get<0>(this->bodies[j]);
			double mj = std::get<2>(this->bodies[j]);
			Eigen::Vector3d rij = xj - xi;
			double dist = rij.norm() + 1e-9; // Vermeide Division durch Null
			force += this->G * mi * mj / (dist * dist * dist) * rij;
		}

		// Euler-Integration
		Eigen::Vector3d vi = std::get<1>(this->bodies[i]);
		Eigen::Vector3d ai = force / mi;
		Eigen::Vector3d vi_new = vi + ai * this->dt;
		Eigen::Vector3d xi_new = xi + vi * this->dt;
		//symplektischer Euler, auch nicht stabiler :(
		//Eigen::Vector3d xi_new = xi + vi_new * this->dt;

		newPositions[i] = xi_new;
		newVelocities[i] = vi_new;
	}

	// Update bodies
	for(size_t i = 0; i < this->bodies.size(); ++i)
	{
		this->bodies[i] = std::make_tuple(newPositions[i], newVelocities[i], std::get<2>(this->bodies[i]));
		positions.push_back(newPositions[i]);
	}

	return positions;
}

void Prog1_2Simulation::mouseEvent(QMouseEvent * e)
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
		cameraElevation = std::fmax(std::fmin(cameraElevation, constants::pi<double> - 0.01), 0.01);

		// tell widget to update itself to account for changed position
		this->update();
	}
}

void Prog1_2Simulation::reset(
	std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, double>> const & bodiesParam,
	/*
		[
			Eigen::Vector3d position,
			Eigen::Vector3d velocity,
			double mass
		]
	*/
	double GParam,
	double dtParam
)
{
	this->G = GParam;
	this->dt = dtParam;
	this->bodies = bodiesParam;

	this->lastTimeNS = this->timer.nsecsElapsed();
	this->timer.restart();
	this->lastSimSteps = 0;

	// Clear the previous textures
	for (auto & texture : this->objectTextures)
	{
		GLuint id = texture.id();
		glDeleteTextures(1, &id);
	}
	this->objectTextures.clear();
	this->objectTextures.resize(bodies.size());

	// Iterate over vector length
	for (unsigned long i = 0; i < bodies.size(); ++i)
	{
		// Generate textures for each body with random colors
		GLuint texId = this->objectTextures[i].id();
		glBindTexture(GL_TEXTURE_2D, texId);

		// Random color pixel
		const uint8_t randomPixel[4] = { static_cast<uint8_t>(rand() % 256), static_cast<uint8_t>(rand() % 256), static_cast<uint8_t>(rand() % 256), 255 };

		qDebug() << "Texture ID:" << texId << "Color:" << randomPixel[0] << randomPixel[1] << randomPixel[2];
		// Upload the texture
		glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8_ALPHA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, randomPixel);

		// Set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		if (GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 1);
	}
}
