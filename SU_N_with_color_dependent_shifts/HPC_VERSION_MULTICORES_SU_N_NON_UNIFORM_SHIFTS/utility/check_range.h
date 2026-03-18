#ifndef GHW_MS_CHECK_RANGE_H
#define GHW_MS_CHECK_RANGE_H

#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <cstddef>

inline void check_range() {
	// Define an unsigned long variable
	unsigned long myVar = 12345678901234567890UL;

	// Print the value
	std::cout << "Value of myVar: " << myVar << std::endl;

	// Check and print the limits of unsigned long
	std::cout << "Maximum value of unsigned long: "
	          << std::numeric_limits<unsigned long>::max() << std::endl;
	std::cout << "Minimum value of unsigned long: "
	          << std::numeric_limits<unsigned long>::min() << std::endl;

	// Verify basic operations
	unsigned long long sum = myVar + 1000000000000000000UL;
	std::cout << "Sum of myVar and 1000000000000000000: " << sum << std::endl;

	// check size of vectors
	std::vector<int> vec_int(10);
	std::cout << "Maximum size of int vector: " << vec_int.max_size() << std::endl;
	if (vec_int.max_size() < static_cast<unsigned long long int>(std::pow(10, 14))) {
		throw std::runtime_error("Machine allows too short int vectors, size<10^14");
	}

	std::vector<double> vec_double(10);
	std::cout << "Maximum size of double vector: " << vec_double.max_size() << std::endl;
	if (vec_double.max_size() < static_cast<unsigned long long int>(std::pow(10, 14))) {
		throw std::runtime_error("Machine allows too short double vectors, size<10^14");
	}

	std::cout << "size_t range:" << std::endl;
	std::cout << "Minimum value: " << 0 << std::endl;
	std::cout << "Maximum value: " << std::numeric_limits<size_t>::max() << std::endl;
	if (std::numeric_limits<size_t>::max() < static_cast<unsigned long long int>(std::pow(10, 14))) {
		throw std::runtime_error("size_t for vector.size() has too small range, <10^14.");
	}

	std::cout << "Range of unsigned long:" << std::endl;
	std::cout << "Minimum value: " << std::numeric_limits<unsigned long>::min() << std::endl;
	std::cout << "Maximum value: " << std::numeric_limits<unsigned long>::max() << std::endl;
	if (std::numeric_limits<unsigned long>::max() < 17446744073709551615UL) {
		throw std::runtime_error("unsigned long has too small range");
	}
	/*
	unsigned long long int threshold = 100000000000000ULL;
	if (N_warm > threshold || N_steps > threshold) {
		throw std::runtime_error("N_warm or iterations are too large, beyond 10^14.");
	}
	 */
}

#endif //GHW_MS_CHECK_RANGE_H