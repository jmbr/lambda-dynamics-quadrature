#include <iostream>
#include <iomanip>

#include <cstdlib>

#include "ensemble-average.h"
#include "observable.h"

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

int main(int argc, char *argv[]) {
  if (argc < 8) {
    std::cerr << "Usage: " << argv[0] << " EPSILON SIGMA Q1 Q2 B"
              << " TEMPERATURE BOX-LENGTH\n";
    return EXIT_FAILURE;
  }

  std::cout << std::setprecision(8);

  const double epsilon{std::strtod(argv[1], nullptr)};
  const double sigma{std::strtod(argv[2], nullptr)};
  const double q1{std::strtod(argv[3], nullptr)};
  const double q2{std::strtod(argv[4], nullptr)};
  const double b{std::strtod(argv[5], nullptr)};
  const double temperature{std::strtod(argv[6], nullptr)};
  const double box_length{std::strtod(argv[7], nullptr)};

  model::Potential potential{epsilon, sigma, q1, q2, b};
  model::EnsembleAverage ensemble_average{potential, temperature, box_length};

  {
    auto observable = std::make_shared<model::LambdaObservable>();
    double error;
    const double value = ensemble_average.compute(observable, &error);
    std::cout << "<lambda> = " << value << " +/- " << error << "\n";
  }

  {
    auto observable = std::make_shared<model::DistanceObservable>();
    double error;
    const double value = ensemble_average.compute(observable, &error);
    std::cout << "<r> = " << value << " +/- " << error << "\n";
  }

  return EXIT_SUCCESS;
}
