#include "Coordinate.hpp"
#include <cmath>
#include <stdexcept>

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

float Coordinate::angle(const std::size_t idx, const std::size_t jdx, const std::size_t kdx) const
{
	const float dxij = xyz_[0][jdx] - xyz_[0][idx];
	const float dyij = xyz_[1][jdx] - xyz_[1][idx];
	const float dzij = xyz_[2][jdx] - xyz_[2][idx];
	const float drij = distance(jdx, idx);
	const float dxjk = xyz_[0][kdx] - xyz_[0][jdx];
	const float dyjk = xyz_[1][kdx] - xyz_[1][jdx];
	const float dzjk = xyz_[2][kdx] - xyz_[2][jdx];
	const float drjk = distance(kdx, jdx);
	const float cost = (dxij * dxjk + dyij * dyjk + dzij * dzjk) / (drij * drjk);
	return std::acos(cost);
}
