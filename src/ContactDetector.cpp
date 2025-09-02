#include "ContactDetector.hpp"
#include "src/util/utility.hpp"
#include <unordered_set>

ContactDetector::result_type ContactDetector::run(
	const Coordinate& coordinate,
	const float cutoff, const std::size_t ignore_num,
	const std::size_t max_contact, const ContactDetector::ff_type& previous_ff) const
{
	std::mt19937_64 engine(seed_);
	const ContactDetector::ff_type& contact_pair_data
		= detect_contact_pairs(coordinate, cutoff, ignore_num, previous_ff, engine);
	if (contact_pair_data.first.size() != contact_pair_data.second.size())
		throw std::runtime_error(
			"[error] The parameter number of indices must be the same as that of native-distance");
	if (contact_pair_data.first.empty() or (contact_pair_data.first.size() <= max_contact))
		return std::make_tuple(
			contact_pair_data.first,
			std::vector<float>(contact_pair_data.first.size(), theta0_),
			std::vector<float>(contact_pair_data.first.size(), phi0_),
			contact_pair_data.second);

	const std::vector<std::size_t>& shuffled_indices
		= Utility::fisher_yates_random_choice(
			max_contact, 0, contact_pair_data.first.size() - 1, engine);
	ContactDetector::result_type retval;
	for (const auto& index : shuffled_indices)
	{
		std::get<0>(retval).push_back(contact_pair_data.first[index]);
		std::get<1>(retval).push_back(theta0_);
		std::get<2>(retval).push_back(phi0_);
		std::get<3>(retval).push_back(contact_pair_data.second[index]);
	}
	return retval;
}

ContactDetector::ff_type ContactDetector::detect_contact_pairs(
	const Coordinate& coordinate, const float cutoff, const std::size_t ignore_num,
	const ContactDetector::ff_type& previous_ff, std::mt19937_64& engine) const
{
	const std::size_t& natom = coordinate.atom_num();
	if (ignore_num < 2)
		throw std::runtime_error("[error] 'ignore_num' must be 2 or more.");
	else if (natom < ignore_num + 1)
		throw std::runtime_error("[error] 'ignore_num' must be smaller than n_atoms.");

	ContactDetector::indices_type contact_indices;
	ContactDetector::param_type   native_distances;
	std::unordered_set<std::size_t> black_list;
	for (std::size_t ipair = 0; ipair < previous_ff.first.size(); ++ipair)
	{
		const auto& pair = previous_ff.first[ipair];
		const auto& native_distance = previous_ff.second[ipair];
		const std::size_t& idx = pair[0];
		const std::size_t& jdx = pair[2];
		const float& dist = coordinate.distance(idx, jdx);
		if (dist < cutoff)
		{
			contact_indices.push_back({idx, idx + 1, jdx, jdx + 1});
			native_distances.push_back(native_distance);
			black_list.insert(idx);
			black_list.insert(jdx);
		}
	}
	for (std::size_t idx = 0; idx < natom - ignore_num - 1; ++idx)
	{
		if (black_list.find(idx) != black_list.end()) continue;
		std::vector<std::size_t> j_indices;
		std::vector<float> dists;
		for (std::size_t jdx = idx + ignore_num; jdx < natom - 1; ++jdx)
		{
			if (black_list.find(jdx) != black_list.end()) continue;
			const float& dist = coordinate.distance(idx, jdx);
			if (dist < cutoff)
			{
				j_indices.push_back(jdx);
				dists.push_back(dist);
			}
		}
		if (j_indices.size() < 1) continue;
		else if (j_indices.size() == 1)
		{
			const std::size_t& jdx = j_indices[0];
			const float& dist = dists[0];
			contact_indices.push_back({idx, idx + 1, jdx, jdx + 1});
			native_distances.push_back(dist);
			black_list.insert(idx);
			black_list.insert(jdx);
		}
		else
		{
			std::uniform_int_distribution<> uni_dist(0, j_indices.size() - 1);
			const std::size_t& ipair = uni_dist(engine);
			const std::size_t& jdx = j_indices[ipair];
			const float& dist = dists[ipair];
			contact_indices.push_back({idx, idx + 1, jdx, jdx + 1});
			native_distances.push_back(dist);
			black_list.insert(idx);
			black_list.insert(jdx);
		}
	}
	return std::make_pair(contact_indices, native_distances);
}

//ContactDetector::result_type ContactDetector::collect_indices_and_parameters(
//	const Coordinate& coordinate, const ContactDetector::indices_type& contact_pairs) const
//{
//	ContactDetector::param_type theta0s;
//	ContactDetector::param_type phi0s;
//	ContactDetector::param_type r0s;
//	for (const auto& pair : contact_pairs)
//	{
//		const float r = coordinate.distance(pair[0], pair[2]);
//		const float theta = coordinate.dihedral(pair[1], pair[0], pair[2], pair[3]);
//		const float phi   = coordinate.angle(pair[1], pair[0], pair[2])
//						  - coordinate.angle(pair[0], pair[2], pair[3]);
//		r0s.push_back(r);
//		if (theta > -0.5 * pi_ && theta < 0.5 * pi_) theta0s.push_back(0.0f);
//		else theta0s.push_back(pi_);
//		if (phi   > -0.5 * pi_ && phi   < 0.5 * pi_) phi0s.push_back(0.0f);
//		else phi0s.push_back(pi_);
//	}
//	return std::make_tuple(contact_pairs, theta0s, phi0s, r0s);
//}
