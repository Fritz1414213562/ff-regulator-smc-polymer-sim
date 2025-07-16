#include "ContactDetector.hpp"
#include "src/util/utility.hpp"
#include <random>
#include <unordered_set>
#include <algorithm>

ContactDetector::result_type ContactDetector::run_mmc(
	const Coordinate& coordinate,
	const float cutoff, const float chi,
	const ContactDetector::indices_type& previous_pairs) const
{
	std::unordered_set<std::size_t> black_list;
	for (const auto& pair : previous_pairs)
	{
		black_list.insert(pair[0]);
		black_list.insert(pair[2]);
	}

	ContactDetector::indices_type indices_vec;
	ContactDetector::param_type   theta0s;
	ContactDetector::param_type   phi0s;

	std::mt19937_64 engine(seed_);
	std::uniform_real_distribution<> dist(0.0, 1.0);
	// try to delete contact pair
	for (const auto& pair : previous_pairs)
	{
		const float dU = - calculate_contact_energy(coordinate, pair) - chi;
		const float accept_rate = std::min(1.0f, std::exp(-dU));
		if (dist(engine) > accept_rate)
		{
			const float theta = coordinate.dihedral(pair[1], pair[0], pair[2], pair[3]);
			const float phi
				= coordinate.angle(pair[1], pair[0], pair[2])
				- coordinate.angle(pair[0], pair[2], pair[3]);
			indices_vec.push_back(pair);
			if (theta > -0.5*pi_ && theta < 0.5*pi_) theta0s.push_back(0.0f);
			else theta0s.push_back(pi_);
			if (phi > -0.5*pi_ && phi < 0.5*pi_) phi0s.push_back(0.0f);
			else phi0s.push_back(pi_);
		}
	}

	// try to append contact pair
	const std::size_t& natom = coordinate.atom_num();
	for (std::size_t idx = 0; idx < natom - 3; ++idx)
	{
		if (black_list.find(idx) != black_list.end()) continue;
		for (std::size_t jdx = idx + 2; jdx < natom - 1; ++jdx)
		{
			if (black_list.find(jdx) != black_list.end()) continue;
			if (coordinate.distance(idx, jdx) < cutoff) 
			{
				const std::array<std::size_t, 4> newpair = {idx, idx + 1, jdx, jdx + 1};
				const float dU = calculate_contact_energy(coordinate, newpair) + chi;
				const float accept_rate = std::min(1.0f, std::exp(-dU));
				if (dist(engine) <= accept_rate)
				{
					const float theta = coordinate.dihedral(idx + 1, idx, jdx, jdx + 1);
					const float phi
						= coordinate.angle(idx + 1, idx, jdx)
						- coordinate.angle(idx, jdx, jdx + 1);
					if (theta > -0.5*pi_ && theta < 0.5*pi_) theta0s.push_back(0.0f);
					else theta0s.push_back(pi_);
					if (phi > -0.5*pi_ && phi < 0.5*pi_) phi0s.push_back(0.0f);
					else phi0s.push_back(pi_);
					indices_vec.push_back(newpair);
					black_list.insert(idx);
					black_list.insert(jdx);
				}
			}
		}
	}
	return std::make_tuple(indices_vec, theta0s, phi0s);
}

float ContactDetector::calculate_contact_energy(const Coordinate& coordinate,
	const std::array<std::size_t, 4>& indices) const
{
	const std::size_t& idx = indices[0];
	const std::size_t& jdx = indices[1];
	const std::size_t& kdx = indices[2];
	const std::size_t& ldx = indices[3];
	const float r = coordinate.distance(idx, kdx);
	const float phi
		= coordinate.angle(jdx, idx, kdx)
		- coordinate.angle(idx, kdx, ldx);
	const float theta = coordinate.dihedral(jdx, idx, kdx, ldx);
//	const float Ub = (r < r0_) ? 0.0 : bond_k_ * (r - r0_) * (r - r0_);
//	const float Ua = angle_k_ * (1.0 - std::cos(2 * (phi - phi0_)));
//	const float Ud = angle_k_ * (1.0 - std::cos(2 * (theta - theta0_)));
//	return Ua + Ub + Ud;
	const float theta0 = (theta > -0.5*pi_ && theta < 0.5*pi_) ? 0.0f : pi_;
	const float phi0   = (phi   > -0.5*pi_ && phi   < 0.5*pi_) ? 0.0f : pi_;
	const float dr = r - r0_;
	const float fr = std::exp(- dr * dr / (2 * sigma_ * sigma_));
	const float K = pi_ / (2 * angle_k_);
	const float dphi = K * (phi - phi0);
	const float dtheta = K * (theta - theta0);
	const float rect0t = (dtheta > -pi_) && (dtheta < pi_) ? 1.0f : 0.0f;
	const float rect1t = (dtheta > -pi_/2) && (dtheta < pi_/2) ? 1.0f : 0.0f;
	const float rect0p = (dphi > -pi_) && (dphi < pi_) ? 1.0f : 0.0f;
	const float rect1p = (dphi > -pi_/2) && (dphi < pi_/2) ? 1.0f : 0.0f;
	const float g1 = 1.0f - std::cos(dtheta) * std::cos(dtheta);
	const float g2 = 1.0f - std::cos(dphi) * std::cos(dphi);
	const float gtheta = std::max(g1 * rect0t, rect1t);
	const float gphi   = std::max(g2 * rect0p, rect1p);
	return -bond_k_ * fr * gtheta * gphi;
}

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
//			if ((coordinate.distance(idx, jdx) < cutoff) and
//				(std::abs(coordinate.angle(idx + 1, idx, jdx) - phi0_) < phi_cutoff_) and
//				(std::abs(coordinate.angle(idx, jdx, jdx + 1) - phi0_) < phi_cutoff_))
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
