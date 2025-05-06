#version 330 core
#define PI 3.1415926535897932384626433832795
#define TWO_PI 6.283185307179586476925286766559

// inputs from vertex shader, different for each fragment
in vec3 objectPosition;
in vec3 viewPosition;
in vec3 viewNormal;

// uniform variables, same for all fragments
uniform sampler2D colorTexture;

// outputs of the fragment shader, i.e., shaded pixels
out vec4 color;

// convert a cartesian unit normal vector to spherical coordinates
vec2 toSphereCoordinate(vec3 normal)
{
	return vec2(atan(normal.y, normal.x) / TWO_PI, 1 - acos(normal.z) / PI);
}

// compute the screen space derivative in x of a cartesian unit normal vector converted to spherical coordinates
vec2 toSphereCoordinateGradX(vec3 n)
{
	// screen space derivative in x
	vec3 ndx = dFdx(n);

	float xy2 = dot(n.xy, n.xy);
	if(xy2 < 1e-6)
		return vec2(0.0, 0.0);
	return vec2((n.x * ndx.y - n.y * ndx.x) / xy2 / TWO_PI, ndx.z / sqrt(1.0 - n.z * n.z) / PI);
}

// compute the screen space derivative in y of a cartesian unit normal vector converted to spherical coordinates
vec2 toSphereCoordinateGradY(vec3 n)
{
	// screen space derivative in y
	vec3 ndx = dFdy(n);

	float xy2 = dot(n.xy, n.xy);
	if(xy2 < 1e-6)
		return vec2(0.0, 0.0);
	return vec2((n.x * ndx.y - n.y * ndx.x) / xy2 / TWO_PI, ndx.z / sqrt(1.0 - n.z * n.z) / PI);
}

void main()
{
	// get surface reflectance (albedo) from texture in spherical coordinates
	vec3 objectNormal = clamp(normalize(objectPosition), -1.0, 1.0);
	vec3 albedo = textureGrad(colorTexture, toSphereCoordinate(objectNormal), toSphereCoordinateGradX(objectNormal), toSphereCoordinateGradY(objectNormal)).rgb;

	// compute simple lambert lighting with light at camera position and no distance-based falloff
	vec3 n = normalize(viewNormal);
	vec3 l = -normalize(viewPosition);
	vec3 lit = albedo * vec3(max(dot(n, l), 0.0));

	// write output color
	color = vec4(lit, 1.0);
}
