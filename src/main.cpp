#include "Input.hpp"
#include "ContactDetector.hpp"
#include "Coordinate.hpp"
#include "src/io/dcd_parser.hpp"
#include "src/io/forcefield_reader.hpp"
#include "src/io/forcefield_writer.hpp"
#include <vector>
#include <array>
#include <unordered_set>
#include <tuple>
#include <utility>


int main(int argc, char *argv[]) {

	using indices_type = std::vector<std::array<std::size_t, 4>>;
	using param_type   = std::vector<float>;
	using ff_type      = std::pair<indices_type, param_type>;
	using result_type  = std::tuple<indices_type, param_type, param_type, param_type>;

	// parse command-line arguments
	Input input(argc, argv);
	// pick up the last frame
	DCDParser dcd_parser = DCDParser();
	Coordinate         coord = dcd_parser.read(input.trajectory_name(), -1);
	// read black list
	ForceFieldReader reader = ForceFieldReader();
	const ff_type& previous_pairs = reader.read(input.base_ff_name());
	// define contact pairs
	ContactDetector detector = ContactDetector(input.seed());
	result_type contact_data
		= detector.run(coord, input.cutoff(), input.ignore_num(), input.max_contact(), previous_pairs);
	// dump the segment-parallelization parameters to the output file
	ForceFieldWriter writer = ForceFieldWriter();
	writer.dump(input.output_name(), std::get<0>(contact_data), input.bond_k(),
		std::get<3>(contact_data), input.sigma(), input.dihedral_k(),
		std::get<1>(contact_data), std::get<2>(contact_data));

	return 0;
}
