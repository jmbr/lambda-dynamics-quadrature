#include <fstream>
#include <iomanip>
#include <iostream>

#include "potential.h"

#include <catch2/catch.hpp>

TEST_CASE("Potential", "[potential]") {
  using namespace model;

  const double L{1e1};
  const size_t n{200};
  const double dr{L / double(n - 1)};
  const double dl{1.0 / double(n - 1)};

  Potential pot;

  std::ofstream out{"potential.dat"};

  out << std::setprecision(8) << std::fixed;

  for (size_t i = 0; i < n; ++i) {
    const double r = i * dr;
    for (size_t j = 0; j < n; ++j) {
      const double lambda = j * dl;
      out << r << "\t" << lambda << "\t" << pot(r, lambda) << "\n";
    }
    out << "\n";
  }
}
