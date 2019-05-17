#include <limits>

#include <cmath>

#include "potential.h"

namespace model {

constexpr double machine_epsilon{std::numeric_limits<double>::epsilon()};

Potential::Potential() : epsilon(1.0), sigma(1.0), q1(-1.0), q2(1.0), b(0.0) {}

Potential::Potential(double epsilon_, double sigma_, double q1_, double q2_,
                     double b_)
    : epsilon(epsilon_), sigma(sigma_), q1(q1_), q2(q2_), b(b_) {}

double Potential::operator()(double r, double lambda) const {
  if (r < machine_epsilon) {
    return INFINITY;
  }

  const double sr2 = sigma / r * sigma / r;
  const double sr6 = sr2 * sr2 * sr2;
  const double sr12 = sr6 * sr6;
  return 4.0 * epsilon * (sr12 - sr6) + lambda * q1 * q2 / r -
         4.0 * b * (lambda - 0.5) * (lambda - 0.5);
}

} // namespace model
