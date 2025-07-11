#ifndef MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
#define MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP

#include <fstream>
#include <iostream>
#include <toml.hpp>
#include <string>
#include <unordered_set>

class ForceFieldReader
{
	public:
		ForceFieldReader() = default;
		~ForceFieldReader() = default;

		std::unordered_set<std::size_t> read(const std::string& filename) const;
};

#endif // MONTE_CARLO_FF_REGULATOR_FORCE_FIELD_READER_HPP
