#include "ExampleRenderer.hpp"

#include <QDebug>
#include <QImage>
#include <QMouseEvent>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <array>

#include <cmath>

#include "helpers.hpp"

ExampleRenderer::ExampleRenderer(QObject * parent)
	: OpenGLRenderer{parent}
{
	{
		// create/bind a vertex array object to store bindings to array and element buffers
		glBindVertexArray(this->icosphereVAO.id());

		// build the required vertex and index arrays
		std::vector<float> vertices(std::begin(icosahedronVertices), std::end(icosahedronVertices));
		std::vector<unsigned> indices(std::begin(icosahedronIndices), std::end(icosahedronIndices));
		for(auto k = 4; k--;)
			subdivideIcosphere(vertices, indices);

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

	{
		// enable seamless cube map filtering
		glEnable(GL_TEXTURE_CUBE_MAP_SEAMLESS);

		// create/bind cube map texture object
		glBindTexture(GL_TEXTURE_CUBE_MAP, this->starsCubeMap.id());

		// load and upload texture images to the individual faces of the cube map
		for (auto nameAndTarget : std::array<std::pair<char const *, GLenum>, 6>{{
			{":/textures/stars_px.jpg", GL_TEXTURE_CUBE_MAP_POSITIVE_X},
			{":/textures/stars_py.jpg", GL_TEXTURE_CUBE_MAP_POSITIVE_Y},
			{":/textures/stars_pz.jpg", GL_TEXTURE_CUBE_MAP_POSITIVE_Z},
			{":/textures/stars_nx.jpg", GL_TEXTURE_CUBE_MAP_NEGATIVE_X},
			{":/textures/stars_ny.jpg", GL_TEXTURE_CUBE_MAP_NEGATIVE_Y},
			{":/textures/stars_nz.jpg", GL_TEXTURE_CUBE_MAP_NEGATIVE_Z}
		}})
		{
			auto img = QImage(nameAndTarget.first).convertToFormat(QImage::Format_RGBA8888)
#if QT_VERSION < QT_VERSION_CHECK(6,9,0)
				.mirrored();
#else
				.flipped(Qt::Vertical);
#endif
			glTexImage2D(nameAndTarget.second, 0, GL_SRGB8_ALPHA8, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV, img.constBits());
		}

		// set texture parameters
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		if(GLAD_GL_EXT_texture_filter_anisotropic)
			glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16);

		// generate mipmaps (required for GL_TEXTURE_MIN_FILTER: GL_LINEAR_MIPMAP_LINEAR)
		glGenerateMipmap(GL_TEXTURE_CUBE_MAP);
	}

	// enable depth tests
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.);

	// enable back face culling
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glFrontFace(GL_CCW);

	this->timer.start();
}

void ExampleRenderer::resize(int w, int h)
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

void ExampleRenderer::render()
{
	auto currentTimeNS = this->timer.nsecsElapsed();
	auto deltaTimeNS = currentTimeNS - this->lastTimeNS;
	this->lastTimeNS = currentTimeNS;

	int loc = -1;

	// clear depth buffer
	glClear(GL_DEPTH_BUFFER_BIT);

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
		Eigen::Affine3d moonModelTransform = Eigen::Translation3d{0., 2., 0.} * Eigen::Scaling(.1);

		// rotate moon around global Z axis with a 10-second period
		auto angle = (currentTimeNS % 10000000000 / 1e10) * constants::two_pi<double>;
		moonModelTransform = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()) * moonModelTransform;

		// bind texture to unit 0 (unit 0 already active)
		glBindTexture(GL_TEXTURE_2D, this->moonTexture.id());

		// set modelView matrix (location already retrieved)
		glUniformMatrix4fv(loc, 1, GL_FALSE, (this->viewMatrix * moonModelTransform.matrix()).cast<float>().eval().data());

		// draw the sphere's triangles as indexed elements (using unsigned int indices), VAO already bound
		glDrawElements(GL_TRIANGLES, this->numIcosphereIndices, GL_UNSIGNED_INT, nullptr);
	}

	// set a different shader
	glUseProgram(this->skyboxProgram.id());

	// set which texture unit is bound to colorTexture (different program, so location may be different too!)
	loc = glGetUniformLocation(this->skyboxProgram.id(), "colorTexture");
	glUniform1i(loc, 0);

	// bind cube map texture to unit 0 (unit 0 already active)
	glBindTexture(GL_TEXTURE_CUBE_MAP, this->starsCubeMap.id());

	// set modelView matrix
	loc = glGetUniformLocation(this->skyboxProgram.id(), "modelView");
	glUniformMatrix4fv(loc, 1, GL_FALSE, this->viewMatrix.cast<float>().eval().data());

	// set projection matrix
	loc = glGetUniformLocation(this->skyboxProgram.id(), "projection");
	glUniformMatrix4fv(loc, 1, GL_FALSE, this->projectionMatrix.cast<float>().eval().data());

	// draw the "skybox"'s triangle as a non-indexed array
	glDepthFunc(GL_EQUAL);
	glBindVertexArray(this->skyboxVAO.id());
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glDepthFunc(GL_LESS);

	// continuous animation, so we always want to redraw
	this->update();
}

void ExampleRenderer::mouseEvent(QMouseEvent * e)
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
