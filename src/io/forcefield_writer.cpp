#include "forcefield_writer.hpp"
#include <iomanip>


void ForceFieldWriter::dump(
	const std::string& filename, const ForceFieldWriter::indices_type& indices_vec,
	const float bond_k, const float r0, const float sigma, const float dihedral_k,
	const ForceFieldWriter::param_type& theta0s, const ForceFieldWriter::param_type& phi0s)
{
	if (indices_vec.size() != theta0s.size() || indices_vec.size() != phi0s.size())
		throw std::runtime_error(
		"[error] The size of contact pairs is not consistent with that of parameters");
	std::ofstream ofs(filename, std::ios::out);
	ofs << "[[forcefields.local]]" << std::endl;
	ofs << "interaction = 'BondLength'" << std::endl;
	ofs << "potential   = 'SegmentParallelization'" << std::endl;
	ofs << "topology    = 'bond'" << std::endl;
	ofs << "parameters  = [" << std::endl;
	const std::size_t natom = indices_vec.size();
	for (std::size_t iatom = 0; iatom < natom; ++iatom)
	{
		const auto& contact_pairs = indices_vec[iatom];
		ofs << "{indices = ["    << std::setw(6) << contact_pairs[0] << ","
		                         << std::setw(6) << contact_pairs[1] << ","
		                         << std::setw(6) << contact_pairs[2] << ","
		                         << std::setw(6) << contact_pairs[3] << "], ";
		ofs << "v0 = "           << std::fixed << std::setprecision(6) <<     r0
			<< ", bond_k = "     << std::fixed << std::setprecision(6) << bond_k
			<< ", theta0 = "     << std::fixed << std::setprecision(6) << theta0s[iatom]
			<< ", phi0 = "       << std::fixed << std::setprecision(6) << phi0s[iatom]
			<< ", sigma = "      << std::fixed << std::setprecision(6) << sigma
			<< ", dihedral_k = " << std::fixed << std::setprecision(6) << dihedral_k << "},"
			<< std::endl;
	}
	ofs << "]" << std::endl;
	ofs.close();
}
