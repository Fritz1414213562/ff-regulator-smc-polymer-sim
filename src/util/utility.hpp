#ifndef MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
#define MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
#include <random>

namespace Utility
{
template<typename T>
T read_binary_as(const char *str)
{
	T result;
	std::memcpy(std::addressof(result), str, sizeof(T));
	return result;
}

std::vector<std::size_t> fisher_yates_random_choice(
	const std::size_t size, std::size_t lower, std::size_t upper, std::mt19937& engine);

}


#endif // MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
