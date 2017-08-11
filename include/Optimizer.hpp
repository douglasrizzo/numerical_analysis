/**
 * @author Douglas De Rizzo Meneghetti (douglasrizzom@gmail.com)
 * @date   2017-8-11
 * @brief Numerical optimizer specialized in finding roots and minima of functions
 */

#ifndef NUMERICAL_ANALYSIS_OPTIMIZER_HPP
#define NUMERICAL_ANALYSIS_OPTIMIZER_HPP

#include <functional>
#include <cmath>
#include <iostream>
#include "FunctionUtils.hpp"

//! Numerical optimizer specialized in finding roots and minima of functions
class Optimizer {
 private:
  int iterations = 0;
  double error = 0;

 public:
  //! \return number of iterations until convergence
  int getIterations() const {
    return iterations;
  }

  //! \return error of approximation
  double getError() const {
    return error;
  }

  //! Numerically searches for the root of a function via Newton-Raphson method
  //! \param f a function
  //! \param x the initial point to start the search
  //! \param error minimum tolerance for the search to end
  //! \param max_iters maximum number of iterations
  //! \param learnRate the learning rate of the search
  //! \param verbose whether to print a short summary of the search at every iteration
  //! \return the point at which the function intercepts the x-axis
  double findRoot(const std::function<double(double)> &f,
                  double x,
                  double error,
                  int max_iters,
                  double learnRate, bool verbose) {
    double f_val = f(x);
    int iter = 0;

    while (iter ++ < max_iters and fabs(f_val) > error) {
      double aux = x + learnRate * - f_val / FunctionUtils::derivative(f, x);
      if (aux == x)
        break;
      x = aux;

      f_val = f(x);
      if (verbose) {
        cout << "Iteration " << iter << ": x = " << x << ", f(x) = " << f_val << '\n';
      }
    }
    this->iterations = iter;
    this->error = fabs(f_val);
    return x;
  }

  //! Function minimization procedure via the gradient descent method
  //! \param f a function
  //! \param x the initial point to start the search
  //! \param error minimum tolerance for the search to end
  //! \param max_iters maximum number of iterations
  //! \param learnRate the learning rate of the search
  //! \param verbose whether to print a short summary of the search at every iteration
  //! \return the point at which the function is minimal
  double minimize(const std::function<double(double)> &f,
                  double x,
                  double error,
                  int max_iters,
                  double learnRate,
                  bool verbose) {
    int iter = 0;
    double d = error + 1;

    while (iter ++ < max_iters and fabs(d) > error) {
      double aux = x - learnRate * d;
      if (aux == x)
        break;
      x = aux;
      d = FunctionUtils::derivative(f, x);

      if (verbose) {
        cout << "Iteration " << iter << ": x = " << x << ", f'(x) = " << d << '\n';
      }
    }
    this->iterations = iter;
    this->error = fabs(d);

    return x;
  }

  //! Two-dimentional function minimization procedure via the gradient descent method
  //! \param f a function
  //! \param x the initial x point to start the search
  //! \param y the initial y point to start the search
  //! \param error minimum tolerance for the search to end
  //! \param max_iters maximum number of iterations
  //! \param learnRate the learning rate of the search
  //! \param verbose whether to print a short summary of the search at every iteration
  //! \return a tuple containing the {x, y} points at which the function is minimal
  std::tuple<double, double> minimize(const std::function<double(double, double)> &f,
                                      double x,
                                      double y,
                                      double error,
                                      int max_iters,
                                      double learnRate,
                                      bool verbose) {
    int iter = 0;
    double d = error + 1;

    while (iter ++ < max_iters and fabs(d) > error) {
      double aux = x - learnRate * d;
      double auy = y - learnRate * d;
      if (aux == x and auy == y)
        break;
      x = aux;
      y = auy;
      d = FunctionUtils::derivative(f, x, y);

      if (verbose) {
        cout << "Iteration " << iter << ": x = " << x << ", y = " << y << ", f'(x, y) = " << d << '\n';
      }
    }
    this->iterations = iter;
    this->error = fabs(d);

    return {x, y};
  }
};

#endif //NUMERICAL_ANALYSIS_OPTIMIZER_HPP
