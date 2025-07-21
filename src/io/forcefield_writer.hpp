#ifndef MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP
#define MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP

#include <array>
#include <vector>
#include <fstream>
#include <cmath>

class ForceFieldWriter
{
		using indices_type = std::vector<std::array<std::size_t, 4>>;
		using param_type = std::vector<float>;
	public:
		ForceFieldWriter() = default;
		~ForceFieldWriter() = default;
		void dump(
			const std::string& filename, const indices_type& indices_vec,
			const float bond_k, const param_type& r0s, const float sigma,
			const float dihedral_k, const param_type& theta0s, const param_type& phi0s);

};

#endif // MONTE_CARLO_FF_REGULATOR_FORCEFIELD_WRITER_HPP
