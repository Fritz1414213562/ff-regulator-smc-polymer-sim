#ifndef MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
#define MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP

#include <fstream>
#include <iostream>
#include <toml.hpp>
#include <string>
#include <vector>
#include <array>
#include <utility>

class ForceFieldReader
{
	using indices_type = std::vector<std::array<std::size_t, 4>>;
	using result_type  = std::pair<indices_type, std::vector<float>>;
	public:
		ForceFieldReader() = default;
		~ForceFieldReader() = default;

		result_type read(const std::string& filename) const;
};

#endif // MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
