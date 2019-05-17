#ifndef OBSERVABLE_H
#define OBSERVABLE_H

namespace model {

class Observable {
public:
  virtual ~Observable() = default;

  virtual double operator()(double r, double lambda) const = 0;
};

} // namespace model

#endif
