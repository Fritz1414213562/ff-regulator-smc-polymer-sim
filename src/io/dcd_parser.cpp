#include "dcd_parser.hpp"
#include "src/util/utility.hpp"


std::vector<Coordinate> DCDParser::read(const std::string& filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
		throw std::runtime_error("[error] File or directory not found: " + filename);

	// 1st block
	const std::string& first_block = read_block(ifs);
	const int frame_num = read_frame_num(first_block);
	const bool has_unitcell = read_unitcell_flag(first_block);
	// 2nd block (not used)
	read_block(ifs);
	// 3rd block
	const std::string& third_block = read_block(ifs);
	const int atom_num  = read_atom_num(third_block);

	std::vector<Coordinate> traj;

	for (std::size_t iframe = 0; iframe < frame_num; ++iframe)
	{
		traj.push_back(Coordinate(read_xyz(ifs, atom_num, has_unitcell)));
	}
	ifs.close();

	return traj;
}

Coordinate DCDParser::read(const std::string& filename, int frame)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
		throw std::runtime_error("[error] File or directory not found: " + filename);
	// 1st block
	const std::string& first_block = read_block(ifs);
	const int frame_num = read_frame_num(first_block);
	const bool has_unitcell = read_unitcell_flag(first_block);
	if (frame >= frame_num)
		throw std::runtime_error(
			"[error] Frame number " + std::to_string(frame) + " is out of range");
	// 2nd block (not used)
	read_block(ifs);
	// 3rd block
	const std::string& third_block = read_block(ifs);
	const int atom_num  = read_atom_num(third_block);

	if (frame < 0) frame = frame_num;

	for (std::size_t iframe = 0; iframe < frame - 1; ++iframe) skip_xyz(ifs, has_unitcell);

	Coordinate coord(read_xyz(ifs, atom_num, has_unitcell));
	ifs.close();
	return coord;
}

std::string DCDParser::read_block(std::ifstream& ifs)
{
	std::int32_t block_size;
	std::vector<char> buffer;

	constexpr int int32_size = sizeof(int32_t);
	ifs.read(reinterpret_cast<char*>(&block_size), int32_size);
	buffer.resize(block_size);
	ifs.read(buffer.data(), block_size);

	std::int32_t check_block_size;
	ifs.read(reinterpret_cast<char*>(&check_block_size), int32_size);

	if (block_size != check_block_size)
		throw std::runtime_error(
			"[error] Inconsistent block size:"
				+ std::to_string(block_size)
				+ " != " + std::to_string(check_block_size));
	return std::string(buffer.begin(), buffer.end());
}

void DCDParser::skip_xyz(std::ifstream& ifs, const bool has_unitcell)
{
	if (has_unitcell) read_block(ifs);
	for (std::size_t idim = 0; idim < 3; ++idim) read_block(ifs);
}

std::array<std::vector<float>, 3> DCDParser::read_xyz(
	std::ifstream& ifs, const int atom_num, const bool has_unitcell)
{
	if (has_unitcell) read_block(ifs);
	const std::string& x_block = read_block(ifs);
	const std::string& y_block = read_block(ifs);
	const std::string& z_block = read_block(ifs);
	std::array<std::vector<float>, 3> xyz = {
		read_coordinates(x_block, atom_num),
		read_coordinates(y_block, atom_num),
		read_coordinates(z_block, atom_num)};
	return xyz;
}

int DCDParser::read_frame_num(const std::string& first_block)
{
	return Utility::read_binary_as<int>(&first_block.at(frame_num_index_));
}

int DCDParser::read_atom_num(const std::string& third_block)
{
	return Utility::read_binary_as<int>(&third_block.at(atom_num_index_));
}

std::vector<float> DCDParser::read_coordinates(const std::string& block, const int atom_num) const
{
	std::vector<float> coordinates(atom_num);
	constexpr int float_size = sizeof(float);
	for (std::size_t iatom = 0; iatom < atom_num; ++iatom)
	{
		std::size_t pos_in_block = iatom * float_size;
		coordinates[iatom] = Utility::read_binary_as<float>(&block.at(pos_in_block));
	}
	return coordinates;
}

bool DCDParser::read_unitcell_flag(const std::string& first_block)
{
	return static_cast<bool>(Utility::read_binary_as<int>(&first_block.at(unit_cell_index_)));
}
