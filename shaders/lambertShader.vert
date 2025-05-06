#version 330 core

// input attributes, different for each vertex
in vec3 position;
in vec3 normal;

// uniform variables, same for all vertices
uniform mat4 projection;
uniform mat4 modelView;

// vertex shader outputs
out vec3 objectPosition;
out vec3 viewPosition;
out vec3 viewNormal;

void main()
{
	// forward position in object/model space to fragment shader
	objectPosition = position;
	// compute and forward position in view space to fragment shader
	vec4 vp = modelView * vec4(position, 1.0);
	viewPosition = vp.xyz;

	// compute normal matrix
	mat3 normalMatrix = transpose(inverse(mat3(modelView)));

	// compute and forward view space normal 
	viewNormal = normalize(normalMatrix * normalize(normal));

	// transform and project position
	gl_Position = projection * vp;
}
