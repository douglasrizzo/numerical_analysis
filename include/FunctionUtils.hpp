/**
 * @author Douglas De Rizzo Meneghetti (douglasrizzom@gmail.com)
 * @date   2017-8-8
 * @brief
 */

#ifndef NUM_ANALYSIS_FUNCTIONUTILS_HPP
#define NUM_ANALYSIS_FUNCTIONUTILS_HPP

#include <functional>

class FunctionUtils {
 public:

  static double derivative(const std::function<double(double)> &f, double x, double delta) {
    return (f(x) - f(x - delta)) / delta;
  }

  static double derivative(const std::function<double(double, double)> &f, double x, double y, double delta) {
    return (f(x, y) - f(x - delta, y - delta)) / delta;
  }

  static double findRoot(const std::function<double(double)> &f,
                         double x,
                         double delta,
                         double error,
                         double learnRate) {
    double f_val = f(x);
    int iter = 0;

    while (fabs(f_val) > error or iter < 1000) {
      x = x + learnRate * (- f_val) * 1 / derivative(f, x, delta);
      f_val = f(x);
      iter ++;
    }

    return x;
  }

  static double minimize(const std::function<double(double)> &f,
                         double x,
                         double delta,
                         double error,
                         double learnRate) {
    int iter = 0;
    double d;

    while (fabs((d = derivative(f, x, delta))) > error or iter < 1000) {
      x = x - learnRate * d;
      iter ++;
    }

    return x;
  }

  static std::tuple<double, double> minimize(const std::function<double(double, double)> &f,
                                             double x,
                                             double y,
                                             double delta,
                                             double error,
                                             double learnRate) {
    int iter = 0;
    double d;

    while (fabs((d = derivative(f, x, y, delta))) > error or iter < 1000) {
      x = x - learnRate * d;
      y = y - learnRate * d;
      iter ++;
    }

    return {x, y};
  }
};

#endif
