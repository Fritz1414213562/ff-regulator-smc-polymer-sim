#ifndef MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
#define MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP

#include <fstream>
#include <iostream>
#include <toml.hpp>
#include <string>
#include <vector>
#include <array>

class ForceFieldReader
{
	using indices_type = std::vector<std::array<std::size_t, 4>>;
	public:
		ForceFieldReader() = default;
		~ForceFieldReader() = default;

		indices_type read(const std::string& filename) const;
};

#endif // MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
