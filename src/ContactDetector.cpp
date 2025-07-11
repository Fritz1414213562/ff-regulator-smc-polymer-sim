#include "ContactDetector.hpp"
#include "src/util/utility.hpp"
#include <random>

ContactDetector::indices_type ContactDetector::run(
	const Coordinate& coordinate, const float cutoff, const std::size_t max_contact) const
{
	const ContactDetector::indices_type& contact_pair_list
		= detect_contact_pairs(coordinate, cutoff);
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
	const Coordinate& coordinate, const float cutoff) const
{
	ContactDetector::indices_type retval;
	const std::size_t& natom = coordinate.atom_num();
	for (std::size_t idx = 0; idx < natom - 3; ++idx)
	{
		for (std::size_t jdx = idx + 2; jdx < natom - 1; ++jdx)
		{
			if (coordinate.distance(idx, jdx) < cutoff)
				retval.push_back({idx, idx + 1, jdx, jdx + 1});
		}
	}
	return retval;
}
