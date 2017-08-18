/**
 * @author Douglas De Rizzo Meneghetti (douglasrizzom@gmail.com)
 * @brief Utility functions for numerical functions
 * @date   2017-8-8
 */

#ifndef NUMERICAL_ANALYSIS_FUNCTIONUTILS_HPP
#define NUMERICAL_ANALYSIS_FUNCTIONUTILS_HPP

#include <functional>
#include <limits>
#include <cmath>

using namespace std;

//! Utility functions for numerical functions
class FunctionUtils {
 private:
  //! A shortcut for the root of the machine epsilon, used for calculating a
  //! small but precise value for h in the derivatives
  static double const sqrtMachineEpsilon;

 public:
  //! Numerically approximates the derivative of a single-variable function
  //! \param f a function
  //! \param x the point at which the derivative is to be calculated
  //! \return the derivative of the function
  static double derivative(const function<double(double)> &f, double x) {
    double h = sqrtMachineEpsilon * x;
    return (f(x + h) - f(x)) / h;
  }

  //! Numerically approximates the partial derivative of a two-dimensional function
  //! \param f a function
  //! \param x the first point at which the derivative is to be calculated
  //! \param y the second point at which the derivative is to be calculated
  //! \param which 0 to differentiate x, otherwise y
  //! \return the partial derivative of the function
  static double partialDerivative(const function<double(double, double)> &f, double x,
                                  double y, int which = 0) {
    double h = sqrtMachineEpsilon;
    if (which == 0) {
      h *= x;
      return (f(x, y) - f(x - h, y)) / h;
    }
    h *= y;
    return (f(x, y) - f(x, y - h)) / h;
  }

  static double rectangleRule(const std::function<double(double)> &f, double a, double b) {
    return (b - a) * f((a + b) / 2);
  }

  static double trapezoidRule(const std::function<double(double)> &f, double a, double b) {
    return (b - a) * (f(a) + f(b)) / 2;
  }

  static double simpsonRule(const std::function<double(double)> &f, double a, double b) {
    return (b - a) * (f(a) + 4 * f((a + b) / 2) + f(b)) / 6;
  }
};

const double FunctionUtils::sqrtMachineEpsilon = sqrt(numeric_limits<double>::epsilon());
#endif
