#include "Input.hpp"
#include "ContactDetector.hpp"
#include "Coordinate.hpp"
#include "src/io/dcd_parser.hpp"
#include "src/io/forcefield_reader.hpp"
#include "src/io/forcefield_writer.hpp"
#include <vector>
#include <array>
#include <unordered_set>


int main(int argc, char *argv[]) {

	using indices_type = std::vector<std::array<std::size_t, 4>>;

	// parse command-line arguments
	Input input(argc, argv);
	// pick up the last frame
	DCDParser dcd_parser = DCDParser();
	Coordinate         coord = dcd_parser.read(input.trajectory_name(), -1);
	// read black list
	ForceFieldReader reader = ForceFieldReader();
	const std::vector<std::array<std::size_t, 4>>& previous_pairs
		= reader.read(input.base_ff_name());
	// define contact pairs
	ContactDetector detector = ContactDetector(input.seed());
	indices_type indices_vec
		= detector.run(coord, input.cutoff(), input.max_contact(), previous_pairs);
	// dump the segment-parallelization parameters to the output file
	ForceFieldWriter writer = ForceFieldWriter();
	writer.dump(input.output_name(),
		indices_vec, input.bond_k(), input.r0(), input.dihedral_k(), input.theta0());

	return 0;
}
