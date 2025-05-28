#pragma once

#include <OpenGLObjects.h>

#include <cstdint>

#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Core>

// helper array with the corner vertices of an icosahedron
constexpr inline float icosahedronVertices[] = {
	0.000000f, -1.000000f, 0.000000f,
	0.723600f, -0.447214f, 0.525720f,
	-0.276386f, -0.447214f, 0.850640f,
	-0.894424f, -0.447214f, 0.000000f,
	-0.276386f, -0.447214f, -0.850640f,
	0.723600f, -0.447214f, -0.525720f,
	0.276386f, 0.447214f, 0.850640f,
	-0.723600f, 0.447214f, 0.525720f,
	-0.723600f, 0.447214f, -0.525720f,
	0.276386f, 0.447214f, -0.850640f,
	0.894424f, 0.447214f, 0.000000f,
	0.000000f, 1.000000f, 0.000000f
};

// helper array of indices which form the triangular faces of an icosahedron
constexpr inline GLubyte icosahedronIndices[] = {
	0, 1, 2,
	1, 0, 5,
	0, 2, 3,
	0, 3, 4,
	0, 4, 5,
	1, 5, 10,
	2, 1, 6,
	3, 2, 7,
	4, 3, 8,
	5, 4, 9,
	1, 10, 6,
	2, 6, 7,
	3, 7, 8,
	4, 8, 9,
	5, 9, 10,
	6, 10, 11,
	7, 6, 11,
	8, 7, 11,
	9, 8, 11,
	10, 9, 11
};

namespace std
{
	// specialization of std::hash for std::pair so that std::unordered_map<std::pair<A,B>, T> can be used
	template<typename A, typename B>
	struct hash<pair<A, B>>
	{
		using argument_type = pair<A, B>;
		using result_type = size_t;

		result_type operator()(argument_type const & p) const noexcept
		{
			return hash<decay_t<A>>{}(p.first) ^ (hash<decay_t<B>>{}(p.second) << 1);
		}
	};
}

// helper function to subdivide a triangular mesh and project it onto the unit sphere
void subdivideIcosphere(std::vector<float> & vertices, std::vector<unsigned> & indices);

// helper function to compute a perspective projection matrix with an infinite far range
Eigen::Matrix4d calculateInfinitePerspective(double minimumFieldOfView, double aspectRatio, double zNear);

// helper function to compute a lookat-matrix (camera at eye, target at center, global up direction; up must be linearly independent of (center - eye))
Eigen::Matrix4d calculateLookAtMatrix(Eigen::Vector3d eye, Eigen::Vector3d center, Eigen::Vector3d up);

// helper function to load a Qt resource as an array of char (bytes)
std::vector<char> loadResource(char const * path);

std::string vectorToString(const Eigen::Vector3d& vec);
