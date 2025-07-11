#include "Coordinate.hpp"
#include <cmath>

Coordinate::Coordinate(const std::array<std::vector<float>, 3>& xyz) : xyz_(xyz)
{
	if (!(xyz[0].size() == xyz[1].size() and xyz[0].size() == xyz[2].size()))
		throw std::runtime_error(
			"[error] The vector sizes of x, y, and z in xyz are not consistent");
	atom_num_ = xyz[0].size();
}

float Coordinate::distance(const std::size_t idx, const std::size_t jdx) const
{
	const float dx = xyz_[0][idx] - xyz_[0][jdx];
	const float dy = xyz_[1][idx] - xyz_[1][jdx];
	const float dz = xyz_[2][idx] - xyz_[2][jdx];
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}
