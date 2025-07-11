#include "ContactDetector.hpp"
#include "src/util/utility.hpp"
#include <random>
#include <unordered_set>

ContactDetector::indices_type ContactDetector::run(
	const Coordinate& coordinate,
	const float cutoff, const std::size_t max_contact,
	const ContactDetector::indices_type& previous_pairs) const
{
	const ContactDetector::indices_type& contact_pair_list
		= detect_contact_pairs(coordinate, cutoff, previous_pairs);
	if (contact_pair_list.empty() or (contact_pair_list.size() <= max_contact))
		return contact_pair_list;

	std::mt19937 engine(seed_);
	const std::vector<std::size_t>& shuffled_indices
		= Utility::fisher_yates_random_choice(
			max_contact, 0, contact_pair_list.size() - 1, engine);
	ContactDetector::indices_type retval;
	for (const auto& index : shuffled_indices)
	{
		retval.push_back(contact_pair_list[index]);
	}
	return retval;
}

ContactDetector::indices_type ContactDetector::detect_contact_pairs(
	const Coordinate& coordinate, const float cutoff,
	const ContactDetector::indices_type& previous_pairs) const
{
	ContactDetector::indices_type retval(previous_pairs);
	std::unordered_set<std::size_t> black_list;
	for (const auto& pair : previous_pairs)
	{
		black_list.insert(pair[0]);
		black_list.insert(pair[2]);
	}
	const std::size_t& natom = coordinate.atom_num();
	for (std::size_t idx = 0; idx < natom - 3; ++idx)
	{
		if (black_list.find(idx) != black_list.end()) continue;
		for (std::size_t jdx = idx + 2; jdx < natom - 1; ++jdx)
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
