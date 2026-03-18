#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <set>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <iomanip>

// Define Real as double for all computations
using Real = double;
using Complex = std::complex<Real>;

// Real number vectors
using Real1DVector = std::vector<Real>;
using Real2DVector = std::vector<Real1DVector>;
using Real3DVector = std::vector<Real2DVector>;
using Real4DVector = std::vector<Real3DVector>;
using Real5DVector = std::vector<Real4DVector>;

// Complex number vectors
using Complex1DVector = std::vector<Complex>;
using Complex2DVector = std::vector<Complex1DVector>;
using Complex3DVector = std::vector<Complex2DVector>;
using Complex4DVector = std::vector<Complex3DVector>;
using Complex5DVector = std::vector<Complex4DVector>;

// Integer vectors
using Int1DVector = std::vector<int>;
using Int2DVector = std::vector<std::vector<int>>;
using Int3DVector = std::vector<std::vector<std::vector<int>>>;
using Int4DVector = std::vector<std::vector<std::vector<std::vector<int>>>>;
using Int5DVector = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>;

#endif // TYPES_H