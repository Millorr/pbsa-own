#include "helpers.hpp"

#include <cmath>

#include <unordered_map>

#include <QFile>

#include <Eigen/Geometry>

void subdivideIcosphere(std::vector<float> & vertices, std::vector<unsigned> & indices)
{
	std::unordered_map<std::pair<unsigned, unsigned>, unsigned> newVertexLookup;
	// reserve number of vertices (= vertices.size() / 3) times 3, as total is multiplied by ~4 per subdivision
	newVertexLookup.reserve(vertices.size());

	auto midpointForEdge = [&] (unsigned first, unsigned second) {
		if(first > second)
			std::swap(first, second);
		auto inserted = newVertexLookup.insert({{first, second}, static_cast<unsigned>(vertices.size() / 3)});
		if(inserted.second)
		{
			Eigen::Map<Eigen::Vector3f> e0{vertices.data() + 3 * first};
			Eigen::Map<Eigen::Vector3f> e1{vertices.data() + 3 * second};
			auto newVertex = (e0 + e1).normalized();
			vertices.insert(std::end(vertices), newVertex.data(), newVertex.data() + 3);
		}
		return inserted.first->second;
	};

	std::vector<unsigned> newIndices;
	newIndices.reserve(4 * indices.size());

	for(std::size_t i = 0; i < indices.size(); i += 3)
	{
		unsigned midpoints[3];
		for(int e = 0; e < 3; ++e)
			midpoints[e] = midpointForEdge(indices[i + e], indices[i + (e + 1) % 3]);
		for(int e = 0; e < 3; ++e)
		{
			newIndices.emplace_back(indices[i + e]);
			newIndices.emplace_back(midpoints[e]);
			newIndices.emplace_back(midpoints[(e + 2) % 3]);
		}
		newIndices.insert(std::end(newIndices), std::begin(midpoints), std::end(midpoints));
	}

	indices.swap(newIndices);
}

Eigen::Matrix4d calculateInfinitePerspective(double minimumFieldOfView, double aspectRatio, double zNear)
{
	// linear field of view factor
	auto range = std::tan(minimumFieldOfView / 2);

	// scale factor for left/right depending on field of view and aspect ratio
	auto right = aspectRatio >= 1.0 ? range * aspectRatio : range;

	// scale factor for top/bottom depending on field of view and aspect ratio
	auto top = aspectRatio >= 1.0 ? range : range / aspectRatio;

	Eigen::Matrix4d P;
	P <<
		1 / right, 0, 0, 0,
		0, 1 / top, 0, 0,
		0, 0, 0, -2 * zNear,
		0, 0, -1, 0;
	return P;
}

Eigen::Matrix4d calculateLookAtMatrix(Eigen::Vector3d eye, Eigen::Vector3d center, Eigen::Vector3d up)
{
	// compute forward direction
	Eigen::RowVector3d f = (eye - center).normalized();
	// compute orthogonal sideways/right direction
	Eigen::RowVector3d s = up.cross(f).normalized();
	// compute orthogonal final up direction
	Eigen::RowVector3d u = f.cross(s);

	// s/u/f define a rotation matrix, inverse translation by rotated eye position, unit row in w for affine transformation
	Eigen::Matrix4d M;
	M <<
		s, -s.dot(eye),
		u, -u.dot(eye),
		f, -f.dot(eye),
		Eigen::RowVector4d::UnitW();
	return M;
}

std::vector<char> loadResource(char const * path)
{
	QFile f(QString(":/") + QString::fromUtf8(path));
	auto opened = f.open(QIODevice::ReadOnly);
	assert(opened); (void)opened;
	auto size = f.size();
	std::vector<char> buf(size);
	auto read = f.read(buf.data(), size);
	assert(read == size); (void)read;
	return buf;
}
