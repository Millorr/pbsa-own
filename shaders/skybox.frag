#version 330 core

// inputs from vertex shader, different for each fragment
in vec3 viewDirection;

// uniform variables, same for all fragments
uniform samplerCube colorTexture;

// outputs of the fragment shader, i.e., shaded pixels
out vec4 color;

void main()
{
	// perform lookup in cube map and output color
	color = texture(colorTexture, viewDirection);
}
