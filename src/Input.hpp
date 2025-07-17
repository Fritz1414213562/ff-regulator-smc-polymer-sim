#ifndef MONTE_CARLO_FF_REGULATOR_INPUT_HPP
#define MONTE_CARLO_FF_REGULATOR_INPUT_HPP

#include <boost/program_options.hpp>
#include <iostream>
#include <string>

class Input {
	public:
		Input(int argc, char *argv[]);
		std::string trajectory_name() const { return trajectory_name_; }
		std::string base_ff_name()    const { return base_ff_name_; }
		std::string output_name()     const { return output_name_; }
		float       bond_k()          const { return bond_k_; }
		float       dihedral_k()      const { return dihedral_k_; }
		float       r0()              const { return r0_; }
		float       sigma()           const { return sigma_; }
		float       cutoff()          const { return cutoff_; }
		int         seed()            const { return seed_; }
		std::size_t max_contact()     const { return max_contact_; }
	private:
		std::string trajectory_name_;
		std::string base_ff_name_;
		std::string output_name_;
		float       bond_k_;
		float       dihedral_k_;
		float       r0_;
		float       sigma_;
		float       cutoff_;
		int         seed_;
		std::size_t max_contact_;
};

#endif // MONTE_CARLO_FF_REGULATOR_INPUT_HPP
