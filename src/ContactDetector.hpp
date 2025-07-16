#ifndef MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
#define MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP

#include "Coordinate.hpp"
#include <vector>
#include <array>
#include <tuple>
#include <cmath>

class ContactDetector
{
	using indices_type = std::vector<std::array<std::size_t, 4>>;
	using param_type = std::vector<float>;
	using result_type = std::tuple<indices_type, param_type, param_type>;

	public:
		ContactDetector(const int seed,
			const float bond_k, const float angle_k, const float r0, const float sigma)
			: seed_(seed), bond_k_(bond_k), angle_k_(angle_k), r0_(r0), sigma_(sigma) {}
		ContactDetector(const int seed) : seed_(seed) {}
		~ContactDetector() = default;
		result_type run_mmc(const Coordinate& coordinate,
			const float cutoff, const float chi,
			const indices_type& previous_pairs) const;
		indices_type run(const Coordinate& coordinate,
			const float cutoff, const std::size_t max_contact,
			const indices_type& previous_pairs) const;

	private:
		indices_type detect_contact_pairs(const Coordinate& coordinate,
			const float cutoff, const indices_type& previous_pairs) const;
		float calculate_contact_energy(const Coordinate& coordinate,
			const std::array<std::size_t, 4>& indices) const;

	private:
		int seed_;
		float bond_k_;
		float angle_k_;
		float r0_;
		float phi0_;
		float theta0_;
		float sigma_;
		float pi_ = std::acos(-1.0);
};


#endif // MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
