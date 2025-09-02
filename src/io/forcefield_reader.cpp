#include "forcefield_reader.hpp"
#include "src/util/utility.hpp"

ForceFieldReader::result_type ForceFieldReader::read(const std::string& filename) const
{
	const toml::value& data = toml::parse(filename);
	const auto& ff = toml::find(data, "forcefields").at(0);
	ForceFieldReader::indices_type indices_set;
	std::vector<float> v0s;
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
					indices_set.push_back(indices);
					const auto& v0 = Utility::find_parameter<float>(param, env, "v0");
					v0s.push_back(v0);
				}
			}
		}
	}
	return std::make_pair(indices_set, v0s);
}
