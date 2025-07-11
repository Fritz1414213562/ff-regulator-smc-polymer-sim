#ifndef MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
#define MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP

#include "Coordinate.hpp"
#include <vector>
#include <array>
#include <unordered_set>

class ContactDetector
{
	using indices_type = std::vector<std::array<std::size_t, 4>>;

	public:
		ContactDetector(const int seed) : seed_(seed) {}
		~ContactDetector() = default;
		indices_type run(const Coordinate& coordinate,
			const float cutoff, const std::size_t max_contact
			const std::unordered_set<std::size_t>& black_list) const;

	private:
		indices_type detect_contact_pairs(const Coordinate& coordinate,
			const std::unordered_set<std::size_t>& black_list, const float cutoff,) const;

	private:
		int seed_;
};


#endif // MONTE_CARLO_FF_REGULATOR_CONTACT_DETECTOR_HPP
