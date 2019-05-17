#include <iostream>
#include <iomanip>

#include <cmath>

#include <boost/math/quadrature/tanh_sinh.hpp>

#include <catch2/catch.hpp>


static inline double f(double x, double y __attribute__((unused))) {
    // return std::exp(-0.5 * (x - 0.5) * (x - 0.5) / sigma2)
    //         / std::sqrt(2.0 * M_PI) ;
    const double k{1e1};
    return 0.5 * k * (x - 0.5) * (x - 0.5) / std::sqrt(x * (1.0 - x));
}


TEST_CASE("One-dimensional quadrature", "[1d-quadrature]") {
  SECTION("One-dimensional quadrature") {
      using quadrature = boost::math::quadrature::tanh_sinh<double>;
      quadrature integrator;

      const double val{integrator.integrate([](double x) { return f(x, 0.0); }, 0.0, 1.0)};
      REQUIRE(val == Approx(1.9634953938));
  }
}
