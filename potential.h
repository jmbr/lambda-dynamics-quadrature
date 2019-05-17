#ifndef POTENTIAL_H
#define POTENTIAL_H

namespace model {

class Potential {
public:
  Potential();
  Potential(double epsilon, double sigma, double q1, double q2, double b);

  double operator()(double r, double lambda) const;

private:
  double epsilon;
  double sigma;
  double q1, q2;
  double b;
};

} // namespace model

#endif
