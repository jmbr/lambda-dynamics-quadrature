#include <iostream>
#include <memory>

#include <catch2/catch.hpp>

#include "ensemble-average.h"
#include "observable.h"
#include "potential.h"

namespace model {

class DistanceObservable : public Observable {
public:
  double operator()(double r, double lambda __attribute__((unused))) const {
    return r;
  }
};

class LambdaObservable : public Observable {
public:
  double operator()(double r __attribute__((unused)), double lambda) const {
    return lambda;
  }
};

} // namespace model

TEST_CASE("EnsembleAverage", "[ensemble-average]") {
  using namespace model;

  SECTION("Lambda") {
    auto observable = std::make_shared<LambdaObservable>();
    Potential potential;
    EnsembleAverage ensemble_average{potential};
    double error;
    const double value = ensemble_average.compute(observable, &error);
    std::cout << "<lambda> = " << value << " +/- " << error << "\n";
  }

  SECTION("Radius") {
    auto observable = std::make_shared<DistanceObservable>();
    Potential potential;
    EnsembleAverage ensemble_average{potential};
    double error;
    const double value = ensemble_average.compute(observable, &error);
    std::cout << "<r> = " << value << " +/- " << error << "\n";
  }
}
