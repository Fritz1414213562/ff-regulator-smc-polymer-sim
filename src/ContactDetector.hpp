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
	using param_type   = std::vector<float>;
	using result_type  = std::tuple<indices_type, param_type, param_type, param_type>;

	public:
		ContactDetector(const int seed) : seed_(seed) {}
		~ContactDetector() = default;
		result_type run(const Coordinate& coordinate,
			const float cutoff, const std::size_t ignore_num,
			const std::size_t max_contact, const indices_type& previous_pairs) const;

	private:
		indices_type detect_contact_pairs(const Coordinate& coordinate,
			const float cutoff, const std::size_t ignore_num,
			const indices_type& previous_pairs) const;
		result_type  collect_indices_and_parameters(
			const Coordinate& coordinate, const indices_type& contact_pairs) const;

	private:
		int seed_;
		float pi_ = std::acos(-1.0);
};


#endif // MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
