#include "forcefield_writer.hpp"
#include <iomanip>


void ForceFieldWriter::dump(
	const std::string& filename, const ForceFieldWriter::indices_type& indices_vec,
	const float bond_k, const float r0,	const float dihedral_k, const float theta0)
{
	std::ofstream ofs(filename, std::ios::app);
	ofs << "[[forcefields.local]]" << std::endl;
	ofs << "interaction = 'BondLength'" << std::endl;
	ofs << "potential   = 'SegmentParallelization'" << std::endl;
	ofs << "topology    = 'bond'" << std::endl;
	ofs << "parameters  = [" << std::endl;
	for (const auto& contact_pairs : indices_vec)
	{
		ofs << "{indices = ["    << std::setw(6) << contact_pairs[0] << ","
		                         << std::setw(6) << contact_pairs[1] << ","
		                         << std::setw(6) << contact_pairs[2] << ","
		                         << std::setw(6) << contact_pairs[3] << "], ";
		ofs << "v0 = "           << std::fixed << std::setprecision(3) <<     r0
			<< ", bond_k = "     << std::fixed << std::setprecision(3) << bond_k
			<< ", theta0 = "     << std::fixed << std::setprecision(3) << theta0
			<< ", dihedral_k = " << std::fixed << std::setprecision(3) << dihedral_k << "},"
			<< std::endl;
	}
	ofs << "]" << std::endl;
	ofs.close();
}
