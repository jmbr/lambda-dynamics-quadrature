#ifndef MODEL_H
#define MODEL_H

#include <memory>

#include "potential.h"

namespace model {

constexpr double DEFAULT_TEMPERATURE{1e0};
constexpr double DEFAULT_BOX_LENGTH{1e1};
constexpr size_t DEFAULT_MAX_EVAL{1000000};
constexpr double DEFAULT_ABS_ERROR{1e-9};
constexpr double DEFAULT_REL_ERROR{1e-9};

class Observable;

class EnsembleAverage {
public:
  EnsembleAverage(const Potential &potential,
                  double temperature = DEFAULT_TEMPERATURE,
                  double box_length = DEFAULT_BOX_LENGTH);

  double compute(const std::shared_ptr<Observable> observable,
                 double *estimated_error = nullptr,
                 size_t max_eval = DEFAULT_MAX_EVAL,
                 double abs_error = DEFAULT_ABS_ERROR,
                 double rel_error = DEFAULT_REL_ERROR);

public:
  const Potential potential;
  const double temperature;
  const double box_length;
  double partition_function, partition_function_error;
  std::shared_ptr<Observable> observable;
};

} // namespace model

#endif
