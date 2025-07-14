#ifndef MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
#define MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP

#include "Coordinate.hpp"
#include <vector>
#include <array>
#include <cmath>

class ContactDetector
{
	using indices_type = std::vector<std::array<std::size_t, 4>>;

	public:
		ContactDetector(const int seed) : seed_(seed) {}
		~ContactDetector() = default;
		indices_type run(const Coordinate& coordinate,
			const float cutoff, const std::size_t max_contact,
			const indices_type& previous_pairs) const;

	private:
		indices_type detect_contact_pairs(const Coordinate& coordinate,
			const float cutoff, const indices_type& previous_pairs) const;

	private:
		int seed_;
		float phi0_ = std::acos(-1.0) / 2.0;
		float phi_cutoff_ = std::acos(-1.0) / 18.0;
};


#endif // MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
