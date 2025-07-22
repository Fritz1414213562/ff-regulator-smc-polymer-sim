#include "ContactDetector.hpp"
#include "src/util/utility.hpp"
#include <random>
#include <unordered_set>

ContactDetector::result_type ContactDetector::run(
	const Coordinate& coordinate,
	const float cutoff, const std::size_t ignore_num,
	const std::size_t max_contact, const ContactDetector::indices_type& previous_pairs) const
{
	const ContactDetector::indices_type& contact_pair_list
		= detect_contact_pairs(coordinate, cutoff, ignore_num, previous_pairs);
	if (contact_pair_list.empty() or (contact_pair_list.size() <= max_contact))
		return collect_indices_and_parameters(coordinate, contact_pair_list);

	std::mt19937_64 engine(seed_);
	const std::vector<std::size_t>& shuffled_indices
		= Utility::fisher_yates_random_choice(
			max_contact, 0, contact_pair_list.size() - 1, engine);
	ContactDetector::indices_type retval;
	for (const auto& index : shuffled_indices)
	{
		retval.push_back(contact_pair_list[index]);
	}
	return collect_indices_and_parameters(coordinate, retval);
}

ContactDetector::indices_type ContactDetector::detect_contact_pairs(
	const Coordinate& coordinate, const float cutoff,
	const std::size_t ignore_num, const ContactDetector::indices_type& previous_pairs) const
{
	const std::size_t& natom = coordinate.atom_num();
	if (ignore_num < 2)
		throw std::runtime_error("[error] 'ignore_num' must be 2 or more.");
	else if (natom < ignore_num + 1)
		throw std::runtime_error("[error] 'ignore_num' must be smaller than n_atoms.");

	ContactDetector::indices_type retval(previous_pairs);
	std::unordered_set<std::size_t> black_list;
	for (const auto& pair : previous_pairs)
	{
		black_list.insert(pair[0]);
		black_list.insert(pair[2]);
	}
	for (std::size_t idx = 0; idx < natom - ignore_num - 1; ++idx)
	{
		if (black_list.find(idx) != black_list.end()) continue;
		for (std::size_t jdx = idx + ignore_num; jdx < natom - 1; ++jdx)
		{
			if (black_list.find(jdx) != black_list.end()) continue;
			if (coordinate.distance(idx, jdx) < cutoff) 
			{
				retval.push_back({idx, idx + 1, jdx, jdx + 1});
				black_list.insert(idx);
				black_list.insert(jdx);
			}
		}
	}
	return retval;
}

ContactDetector::result_type ContactDetector::collect_indices_and_parameters(
	const Coordinate& coordinate, const ContactDetector::indices_type& contact_pairs) const
{
	ContactDetector::param_type theta0s;
	ContactDetector::param_type phi0s;
	ContactDetector::param_type r0s;
	for (const auto& pair : contact_pairs)
	{
		const float r = coordinate.distance(pair[0], pair[2]);
		const float theta = coordinate.dihedral(pair[1], pair[0], pair[2], pair[3]);
		const float phi   = coordinate.angle(pair[1], pair[0], pair[2])
						  - coordinate.angle(pair[0], pair[2], pair[3]);
		r0s.push_back(r);
		if (theta > -0.5 * pi_ && theta < 0.5 * pi_) theta0s.push_back(0.0f);
		else theta0s.push_back(pi_);
		if (phi   > -0.5 * pi_ && phi   < 0.5 * pi_) phi0s.push_back(0.0f);
		else phi0s.push_back(pi_);
	}
	return std::make_tuple(contact_pairs, theta0s, phi0s, r0s);
}
