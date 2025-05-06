#version 330 core

// corners of a full-screen triangle at the far plane in normalized device coordinates
const vec4 positions[3] = vec4[3](
	vec4(-1.0, -1.0, 1.0, 1.0),
	vec4( 3.0, -1.0, 1.0, 1.0),
	vec4(-1.0,  3.0, 1.0, 1.0)
);

// uniform variables, same for all vertices
uniform mat4 projection;
uniform mat4 modelView;

// vertex shader outputs
out vec3 viewDirection;

void main()
{
	// compute necessary inverse matrices
	mat4 inverseProjection = inverse(projection);
	mat3 inverseNormalMatrix = transpose(mat3(modelView));

	// look up normalized device coordinates from constant using vertex index
	vec4 ndc = positions[gl_VertexID];

	// unproject NDC coordinates to compute a direction
	vec3 unprojected = (inverseProjection * ndc).xyz;
	viewDirection = inverseNormalMatrix * unprojected;

	// projection copied as is
	gl_Position = ndc;
}
