#ifndef MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP
#define MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP

#include <array>
#include <vector>
#include <fstream>

class ForceFieldWriter
{
		using indices_type = std::vector<std::array<std::size_t, 4>>;
	public:
		ForceFieldWriter() = default;
		~ForceFieldWriter() = default;
		void dump(
			const std::string& filename, const indices_type& indices_vec,
			const float bond_k, const float r0,
			const float dihedral_k, const float theta0);
};

#endif // MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP
