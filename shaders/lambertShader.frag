#version 330 core

// inputs from vertex shader, different for each fragment
in vec3 objectPosition;
in vec3 viewPosition;
in vec3 viewNormal;
in vec3 vertexColor;

// uniform variables, same for all fragments
uniform sampler2D colorTexture;

// outputs of the fragment shader, i.e., shaded pixels
out vec4 color;

void main()
{
	// compute simple lambert lighting with light at camera position and no distance-based falloff
	vec3 n = normalize(viewNormal);
	vec3 l = -normalize(viewPosition);
	vec3 c = vertexColor;
	vec3 lit = c * vec3(min(max(dot(n, l), 0.0), 0.5));
	// write output color
	color = vec4(lit, 1.0);
}
