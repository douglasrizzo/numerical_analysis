/**
 * @author Douglas De Rizzo Meneghetti (douglasrizzom@gmail.com)
 * @date   2017-8-8
 * @brief Utility functions for numerical functions
 */

#ifndef NUMERICAL_ANALYSIS_FUNCTIONUTILS_HPP
#define NUMERICAL_ANALYSIS_FUNCTIONUTILS_HPP

#include <functional>
#include <limits>

using namespace std;

//! Utility functions for numerical functions
class FunctionUtils {
 private:
  //! A shortcut for the root of the machine epsilon, used in the calculating a small but precise value for h in the derivatives
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

  //! Numerically approximates the derivative of a two-dimensional function
  //! \param f a function
  //! \param x the first point at which the derivative is to be calculated
  //! \param y the second point at which the derivative is to be calculated
  //! \return the derivative of the function
  static double derivative(const function<double(double, double)> &f, double x, double y) {
    double h = sqrtMachineEpsilon * x;
    return (f(x + h, y + h) - f(x - h, y - h)) / 2 * h;
  }
};

const double FunctionUtils::sqrtMachineEpsilon = sqrt(numeric_limits<double>::epsilon());
#endif
