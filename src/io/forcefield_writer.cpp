#include "forcefield_writer.hpp"
#include <iomanip>


void ForceFieldWriter::dump(
	const std::string& filename, const std::vector<std::array<std::size_t, 4>>& indices_vec,
	const float bond_k, const float r0,	const float dihedral_k, const float sigma,
	const std::vector<float>& phi0, const std::vector<float>& theta0)
{
	if (indices_vec.size() != phi0.size() || indices_vec.size() != theta0.size())
		throw std::runtime_error(
			"[error] The indices size is not consistent with that of theta0 or phi0");
	std::ofstream ofs(filename, std::ios::out);
	ofs << "[[forcefields.local]]" << std::endl;
	ofs << "interaction = 'BondLength'" << std::endl;
	ofs << "potential   = 'SegmentParallelization'" << std::endl;
	ofs << "topology    = 'bond'" << std::endl;
	ofs << "parameters  = [" << std::endl;

	const std::size_t natom = indices_vec.size();
	for (std::size_t iatom = 0; iatom < natom; ++iatom)
	{
		const std::array<std::size_t, 4>& contact_pairs = indices_vec[iatom];
		ofs << "{indices = ["    << std::setw(6) << contact_pairs[0] << ","
		                         << std::setw(6) << contact_pairs[1] << ","
		                         << std::setw(6) << contact_pairs[2] << ","
		                         << std::setw(6) << contact_pairs[3] << "], ";
		ofs << "v0 = "           << std::fixed << std::setprecision(6) <<     r0
			<< ", bond_k = "     << std::fixed << std::setprecision(6) << bond_k
			<< ", theta0 = "     << std::fixed << std::setprecision(6) << theta0[iatom]
			<< ", phi0 = "       << std::fixed << std::setprecision(6) << phi0[iatom]
			<< ", sigma = "      << std::fixed << std::setprecision(6) << sigma
			<< ", dihedral_k = " << std::fixed << std::setprecision(6) << dihedral_k << "},"
			<< std::endl;
	}
	ofs << "]" << std::endl;
	ofs.close();
}
