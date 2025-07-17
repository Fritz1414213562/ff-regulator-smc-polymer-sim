#ifndef MONTE_CARLO_FF_REGULATOR_DCD_PARSER_HPP
#define MONTE_CARLO_FF_REGULATOR_DCD_PARSER_HPP

#include "src/Coordinate.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <array>

class DCDParser {

	public:
		DCDParser()  = default;
		~DCDParser() = default;

		std::vector<Coordinate> read(const std::string& filename);
		Coordinate              read(const std::string& filename, int frame);

	private:
		std::string read_block(std::ifstream& ifs);
		std::array<std::vector<float>, 3> read_xyz(
			std::ifstream& ifs, const int atom_num, const bool has_unitcell);
		void skip_xyz(std::ifstream& ifs, const bool has_unitcell);

		int read_frame_num(const std::string& first_block);
		int read_atom_num(const std::string& third_block);
		bool read_unitcell_flag(const std::string& first_block);
		std::vector<float> read_coordinates(const std::string& block, const int atom_num) const;

	private:
		const int frame_num_index_ = 4;
		const int unit_cell_index_ = 44;
		const int atom_num_index_  = 0;
		const float A2nm = 0.1;

};

#endif // MONTE_CARLO_FF_REGULATOR_DCD_PARSER_HPP
