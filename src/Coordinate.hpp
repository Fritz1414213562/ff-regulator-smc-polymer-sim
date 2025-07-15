#ifndef MONTE_CARLO_FF_REGULATOR_COORDINATE_HPP
#define MONTE_CARLO_FF_REGULATOR_COORDINATE_HPP

#include <vector>
#include <array>

class Coordinate
{
	public:
		Coordinate(const std::array<std::vector<float>, 3>& xyz);
		std::array<std::vector<float>, 3>  xyz() const {return xyz_;}
		std::array<std::vector<float>, 3>& xyz()       {return xyz_;}
		std::size_t                   atom_num() const {return atom_num_;}
		std::size_t&                  atom_num()       {return atom_num_;}
	
		float distance(const std::size_t idx, const std::size_t jdx) const;
		float angle(const std::size_t idx, const std::size_t jdx, const std::size_t kdx) const;
		float dihedral(const std::size_t idx, const std::size_t jdx,
			const std::size_t kdx, const std::size_t ldx) const;

	private:
		std::size_t atom_num_;
		std::array<std::vector<float>, 3> xyz_;
};

#endif // MONTE_CARLO_FF_REGULATOR_COORDINATE_HPP
