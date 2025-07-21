#ifndef MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
#define MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
#include <random>
#include <cstring>
#include <toml.hpp>
#include <toml_fwd.hpp>

namespace Utility
{
template<typename T>
T read_binary_as(const char *str)
{
	T result;
	std::memcpy(std::addressof(result), str, sizeof(T));
	return result;
}

template<typename T>
T find_parameter(const toml::value& params, const toml::value& env,
                 const std::string& name)
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    if(!params.is_table() || !params.contains(name))
    {
        const toml::error_info err = toml::make_error_info(
                "value " + name + " does not exists", params, "in this table");
        throw std::out_of_range(toml::format_error(err));
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_empty())
        {
            const toml::error_info err = toml::make_error_info(
                    "named variable \"" + var + "\" used but no env is defined",
                    params, "used here");
            throw std::out_of_range(toml::format_error(err));
        }
        if(!env.is_table() || !env.contains(var))
        {
            const toml::error_info err = toml::make_error_info(
                    "named variable \"" + var + "\" does not exists",
                    env, "in this table");
            throw std::out_of_range(toml::format_error(err));
        }
        return toml::find<T>(env, var);
    }
    return toml::get<T>(p);
}

template<typename T>
T inner_product(const std::array<T, 3>& lhs, const std::array<T, 3>& rhs)
{
	return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

template<typename T>
std::array<T, 3> cross_product(const std::array<T, 3>& lhs, const std::array<T, 3>& rhs)
{
	return std::array<T, 3>({
		lhs[1] * rhs[2] - lhs[2] * rhs[1],
		lhs[2] * rhs[0] - lhs[0] * rhs[2],
		lhs[0] * rhs[1] - lhs[1] * rhs[0],
	});
}

std::vector<std::size_t> fisher_yates_random_choice(
	const std::size_t size, std::size_t lower, std::size_t upper, std::mt19937_64& engine);

}


#endif // MONTE_CARLO_FF_REGULATOR_UTILITY_HPP
