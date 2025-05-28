#include "Prog1_1Simulation.hpp"

#include "helpers.hpp"

#include <QImage>
#include <QDebug>
#include <QMouseEvent>
#include <QTimer>

#include <Eigen/Dense>
#include <Eigen/Geometry>

Prog1_1Simulation::Prog1_1Simulation()
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
		// load a texture from a Qt resource
		auto img = QImage(":/textures/earth_color.jpg").convertToFormat(QImage::Format_RGBA8888)
#if QT_VERSION < QT_VERSION_CHECK(6,9,0)
			.mirrored();
#else
			.flipped(Qt::Vertical);
#endif

		// create/bind texture object
		glBindTexture(GL_TEXTURE_2D, this->earthTexture.id());
		
		// upload texture data
		glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8_ALPHA8, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV, img.constBits());

		// set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		if(GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16);

		// generate mipmaps (required for GL_TEXTURE_MIN_FILTER: GL_LINEAR_MIPMAP_LINEAR)
		glGenerateMipmap(GL_TEXTURE_2D);
	}

	{
		// load a texture from a Qt resource
		auto img = QImage(":/textures/moon_color.jpg").convertToFormat(QImage::Format_RGBA8888)
#if QT_VERSION < QT_VERSION_CHECK(6,9,0)
			.mirrored();
#else
			.flipped(Qt::Vertical);
#endif

		// create/bind texture object
		glBindTexture(GL_TEXTURE_2D, this->moonTexture.id());

		// upload texture data
		glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB8_ALPHA8, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV, img.constBits());

		// set texture parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		if(GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16);

		// generate mipmaps (required for GL_TEXTURE_MIN_FILTER: GL_LINEAR_MIPMAP_LINEAR)
		glGenerateMipmap(GL_TEXTURE_2D);
	}

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	this->timer.start();
}

Prog1_1Simulation::~Prog1_1Simulation() = default;

Eigen::Vector3d Prog1_1Simulation::explicitEuler()
{
	Eigen::Vector3d a;
	Eigen::Vector3d rt;
	Eigen::Vector3d drt;
	Eigen::Matrix3d A;
	double r_norm3 = this->r.lpNorm<3>();
	a = -this->g / r_norm3 * this->r;
	drt = this->dr + a * this->dt;
	rt = this->r + this->dr * this->dt;
	this->dr = drt;
	this->r = rt;

	return rt;
}

Eigen::Vector3d Prog1_1Simulation::implicitEuler()
{
	// Implicit Euler integration
	Eigen::Vector3d rt = this->r;
	Eigen::Vector3d drt = this->dr;

	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d A;
	Eigen::Matrix3d outer;
	Eigen::Vector3d b;
	Eigen::Vector3d r_next;
	Eigen::Vector3d dr_next;
	constexpr double epsilon = 1e-8;

	double r_norm3 = std::max(rt.lpNorm<3>(), epsilon);
	double r_norm5 = std::max(rt.lpNorm<5>(), epsilon);
	outer = rt * rt.transpose();

	A = I + this->dt * this->dt * this->g * (I / r_norm3 - 3.0 * outer / r_norm5);
	b = rt + this->dt * drt;

	if(std::abs(A.determinant()) < epsilon)
	{
		qDebug() << "Warnung: Matrix A ist fast singulär!";
	}

	rt = A.fullPivLu().solve(b);
	drt = drt - this->dt * this->g / r_norm3 * rt;

	this->r = rt;
	this->dr = drt;
	return rt;
}

Eigen::Vector3d Prog1_1Simulation::verlet()
{
	// Verlet-Integration: x(t+dt) = 2*x(t) - x(t-dt) + a(t)*dt²

	// Berechnung der aktuellen Beschleunigung
	double r_norm3 = this->r.lpNorm<3>();
	Eigen::Vector3d a = -this->g / r_norm3 * this->r;

	Eigen::Vector3d r_next;

	if(this->firstStep)
	{
		// Erster Schritt: Verwende expliziten Euler zur Initialisierung
		// x(t-dt) ≈ x(t) - v(t)*dt
		this->r_prev = this->r - this->dr * this->dt;
		this->firstStep = false;
	}

	// Verlet-Formel: x(t+dt) = 2*x(t) - x(t-dt) + a(t)*dt²
	r_next = 2.0 * this->r - this->r_prev + a * this->dt * this->dt;

	// Update für nächsten Schritt
	this->r_prev = this->r;
	this->r = r_next;

	return r_next;
}

Eigen::Vector3d Prog1_1Simulation::step()
{
	switch(this->integration)
	{
	default:
	case Integration::Explicit:
		// Explicit Euler integration
		return explicitEuler();
		break;
	case Integration::Implicit:
		// Implicit Euler integration
		return implicitEuler();
		break;
	case Integration::Verlet:
		return verlet();
		break;
	}
}

void Prog1_1Simulation::resize(int w, int h)
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

void Prog1_1Simulation::render()
{
	auto currentTimeNS = this->timer.nsecsElapsed();
	auto deltaTimeNS = currentTimeNS - this->lastTimeNS;
	this->lastTimeNS = currentTimeNS;
	qint64 simSteps = (currentTimeNS/this->dt/ 1e9) - this->lastSimSteps;
	this->lastSimSteps += simSteps;

	int loc = -1;

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

	{
		auto tiltAngle = 23.4 * constants::deg_to_rad<double>;
		Eigen::Affine3d earthModelTransform = Eigen::AngleAxisd{tiltAngle, Eigen::Vector3d::UnitY()} * Eigen::Affine3d::Identity();

		// rotate earth around local Z axis with a 3-second period
		auto angle = (currentTimeNS % 3000000000 / 3e9) * constants::two_pi<double>;
		earthModelTransform = earthModelTransform * Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());

		// bind texture to unit 0
		glActiveTexture(GL_TEXTURE0 + 0);
		glBindTexture(GL_TEXTURE_2D, this->earthTexture.id());

		// set modelView matrix
		loc = glGetUniformLocation(this->icosphereProgram.id(), "modelView");
		glUniformMatrix4fv(loc, 1, GL_FALSE, (this->viewMatrix * earthModelTransform.matrix()).cast<float>().eval().data());

		// draw the sphere's triangles as indexed elements (using unsigned int indices), VAO already bound
		glDrawElements(GL_TRIANGLES, this->numIcosphereIndices, GL_UNSIGNED_INT, nullptr);
	}

	{
		Eigen::Affine3d moonModelTransform = Eigen::Scaling(.1) * Eigen::Affine3d::Identity();
		Eigen::Vector3d simPos;
		for(int i = 0; i < simSteps; i++)
		{
			simPos = this->step();
		}
		qDebug() << "FPS: " << 1.0 / (deltaTimeNS / 1.0e09) << " simsteps: " << simSteps << " x: " << simPos.x() << " y: " << simPos.y() << " z: " << simPos.z();

		Eigen::Affine3d simTransform = Eigen::Translation3d{simPos.x(), simPos.y(), simPos.z()} * Eigen::Affine3d::Identity();

		// use simulation data for new position
		moonModelTransform = moonModelTransform * simTransform;

		// bind texture to unit 0 (unit 0 already active)
		glBindTexture(GL_TEXTURE_2D, this->moonTexture.id());

		// set modelView matrix (location already retrieved)
		glUniformMatrix4fv(loc, 1, GL_FALSE, (this->viewMatrix * moonModelTransform.matrix()).cast<float>().eval().data());

		// draw the sphere's triangles as indexed elements (using unsigned int indices), VAO already bound
		glDrawElements(GL_TRIANGLES, this->numIcosphereIndices, GL_UNSIGNED_INT, nullptr);
	}

	// continuous animation, so we always want to redraw
	this->update();
}

void Prog1_1Simulation::mouseEvent(QMouseEvent * e)
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

void Prog1_1Simulation::reset(
	Integration integrationParam,
	Eigen::Vector3d r0Param,
	Eigen::Vector3d dr0Param,
	double gParam,
	double dtParam
)
{

	// Set initial conditions
	this->integration = integrationParam;
	this->r0 = r0Param;
	this->dr0 = dr0Param;
	this->g = gParam;
	this->dt = dtParam;

	// Reset simulation state
	this->r = r0Param;
	this->dr = dr0Param;


	this->lastTimeNS = this->timer.nsecsElapsed();
	this->timer.restart();
	this->lastSimSteps = 0;
	this->firstStep = true;

	qDebug() << "Reset simulation with parameters:"
		<< "integration:" << static_cast<int>(integration)
		<< "r0:" << vectorToString(r0)
		<< "dr0:" << vectorToString(dr0)
		<< "g:" << g
		<< "dt:" << dt;
}
