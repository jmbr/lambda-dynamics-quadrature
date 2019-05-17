#include <algorithm>
#include <iostream>
#include <limits>

#include <cassert>
#include <cmath>
#include <cstdlib>

#include "cubature.h"
#include <boost/math/quadrature/tanh_sinh.hpp>

#include "ensemble-average.h"
#include "observable.h"

namespace model {

using quadrature = boost::math::quadrature::tanh_sinh<double>;

EnsembleAverage::EnsembleAverage(const Potential &potential_,
                                 double temperature_, double box_length_)
    : potential(potential_), temperature(temperature_), box_length(box_length_),
      partition_function(NAN), partition_function_error(INFINITY) {}

// Compute distance taking periodic boundary conditions into account.
static inline double distance(double x, double y, double length) {
  return std::min(
      {std::fabs(x - y), std::fabs(x - y - length), std::fabs(x - y + length)});
}

static inline double jacobian(double lambda) {
  return 1.0 / std::sqrt(lambda * (1.0 - lambda));
}

static inline double integrand_lambda(double lambda, double r,
                                      const EnsembleAverage &avg) {
  return std::exp(-avg.potential(r, lambda) / avg.temperature) *
         jacobian(lambda);
}

static int int_partition_function(unsigned int ndim __attribute__((unused)),
                                  const double *x, void *fdata,
                                  unsigned int fdim __attribute__((unused)),
                                  double *fval) {
  auto &avg = *((EnsembleAverage *)fdata);
  const double r{distance(x[0], x[1], avg.box_length)};

  quadrature integrator;

  *fval = integrator.integrate(
      [&](double lambda) { return integrand_lambda(lambda, r, avg); }, 0.0,
      1.0);

  return 0;
}

static int int_ensemble_average(unsigned int ndim __attribute__((unused)),
                                const double *x, void *fdata,
                                unsigned int fdim __attribute__((unused)),
                                double *fval) {
  auto &avg = *((EnsembleAverage *)fdata);
  const double r{distance(x[0], x[1], avg.box_length)};

  quadrature integrator;

  *fval = integrator.integrate(
      [&](double lambda) {
        return avg.observable->operator()(r, lambda) *
               integrand_lambda(lambda, r, avg) / avg.partition_function;
      },
      0.0, 1.0);

  return 0;
}

double EnsembleAverage::compute(std::shared_ptr<Observable> const observable,
                                double *estimated_error, size_t max_eval,
                                double abs_error, double rel_error) {
  const std::vector<double> xmin{0.0, 0.0};
  const std::vector<double> xmax{box_length, box_length};
  assert(xmin.size() == xmax.size());
  const size_t fdim{1};

  auto cubature = [&](::integrand func) {
    int status;
    double value, error;
    if ((status = hcubature(fdim, func, static_cast<void *>(this), xmin.size(),
                            &xmin[0], &xmax[0], max_eval, abs_error, rel_error,
                            ERROR_INDIVIDUAL, &value, &error)) != 0) {
      throw std::runtime_error{"Error computing integral"};
    }

    return std::make_pair(value, error);
  };

  if (std::isnan(partition_function)) {
    std::tie(partition_function, partition_function_error) =
        cubature(int_partition_function);

    std::cout << "Partition function = " << partition_function << " +/- "
              << partition_function_error << "\n";

    if (partition_function_error / std::fabs(partition_function) > 1e-3) {
      std::cerr << "The value of the partition function, " << partition_function
                << " +/- " << partition_function_error
                << ", might not be accurate enough.\n";
    }
  }

  this->observable = observable;
  double ensemble_average, ensemble_average_error;
  std::tie(ensemble_average, ensemble_average_error) =
      cubature(int_ensemble_average);

  if (estimated_error != nullptr)
    *estimated_error = ensemble_average_error;

  return ensemble_average;
}

} // namespace model
