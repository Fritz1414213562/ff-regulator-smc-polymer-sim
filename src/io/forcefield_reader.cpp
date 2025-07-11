#include "forcefield_reader.hpp"

std::unordered_set<std::size_t> ForceFieldReader::read(const std::string& filename)
{
	const toml::value& data = toml::parse(filename);
	const auto& ff = toml::find(data, "forcefields").at(0);
	std::unordered_set<std::size_t> retval;
	if (ff.contains("local"))
	{
		const auto& locals = toml::find(ff, "local").as_array();
		for (const auto& local_ff : locals)
		{
			const std::string& interaction = toml::find<std::string>(local_ff, "interaction");
			const std::string& potential   = toml::find<std::string>(local_ff, "potential");
			if (interaction == "BondLength" && potential == "SegmentParallelization")
			{
				const auto& params = toml::find<toml::array>(local_ff, "parameters");
				const auto& env = local_ff.contains("env") ? local_ff.at("env") : toml::value{};
				for (const auto& param : params)
				{
					const auto& indices = Utility::find_parameter<std::array<std::size_t, 4>>(
						param, env, "indices");
					retval.insert(indices[0]);
					retval.insert(indices[2]);
				}
			}
		}
	}
	return retval;
}
