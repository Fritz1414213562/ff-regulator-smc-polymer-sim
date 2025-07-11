#include "utility.hpp"
#include <cstring>
#include <memory>
#include <random>
#include <unordered_map>

namespace Utility
{

std::vector<std::size_t> fisher_yates_random_choice(
	const std::size_t size, std::size_t lower, std::size_t upper, std::mt19937& engine)
{
	if (lower > upper) std::swap(lower, upper);
	if (upper - lower + 1 < size) throw std::runtime_error(
		"[error] The range of random choice must be larger than the sample size");

	std::vector<std::size_t> retval;
	retval.reserve(size);

	std::unordered_map<std::size_t, std::size_t> replaced_map;

	for (std::size_t _cycle = 0; _cycle < size; ++_cycle)
	{
		const std::size_t val = std::uniform_int_distribution<std::size_t>(lower, upper)(engine);
		auto itr = replaced_map.find(val);
		std::size_t replaced_val;
		auto replaced_itr = replaced_map.find(upper);

		if (replaced_itr != replaced_map.end()) replaced_val = replaced_itr->second;
		else replaced_val = upper;

		if (itr == replaced_map.end())
		{
			retval.push_back(val);
			if (val != upper) replaced_map.insert(std::make_pair(val, replaced_val));
		}
		else
		{
			retval.push_back(itr->second);
			itr->second = replaced_val;
		}
		--upper;
	}

	return retval;
}

}
