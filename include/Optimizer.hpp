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
  string endReason = "You didn't run any optimization yet!";

 public:
  //! \return number of iterations until convergence
  int getIterations() const {
    return iterations;
  }

  //! \return error of approximation
  double getError() const {
    return error;
  }

  const string &getEndReason() const {
    return endReason;
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
                  double error = 1e-8,
                  int max_iters = 1000,
                  double learnRate = 1, bool verbose = false) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double f_val = f(x);

    while (true) {
      double aux = x + learnRate * - f_val / FunctionUtils::derivative(f, x);
      if (aux == x) {
        this->endReason = "No change in x from previous iteration";
        break;
      }

      iterations ++;
      x = aux;
      f_val = f(x);
      if (verbose) {
        cout << "Iteration " << iterations << ": x = " << x << ", f(x) = " << f_val << '\n';
      }

      if (iterations >= max_iters) {
        this->endReason = "Maximum number of iterations reached";
        break;
      }
      if (fabs(f_val) < error) {
        this->endReason = "Minimum error threshold reached";
        break;
      }
    }
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
                  double error = 1e-8,
                  int max_iters = 1000,
                  double learnRate = 1,
                  bool verbose = false) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double d = error + 1;

    while (true) {
      double aux = x - learnRate * d;
      if (aux == x) {
        this->endReason = "No change in x from previous iteration";
        break;
      }
      iterations ++;
      x = aux;
      d = FunctionUtils::derivative(f, x);

      if (verbose) {
        cout << "Iteration " << iterations << ": x = " << x << ", f'(x) = " << d << '\n';
      }

      if (iterations >= max_iters) {
        this->endReason = "Maximum number of iterations reached";
        break;
      }
      if (fabs(d) < error) {
        this->endReason = "Minimum error threshold reached";
        break;
      }
    }
    this->error = fabs(d);

    return x;
  }

  //! Two-dimensional function minimization procedure via the gradient descent method
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
                                      double error = 1e-8,
                                      int max_iters = 1000,
                                      double learnRate = 1,
                                      bool verbose = false) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double d = error + 1;

    while (true) {
      double aux = x - learnRate * d;
      double auy = y - learnRate * d;
      if (aux == x and y == auy) {
        this->endReason = "No change in x and y from previous iteration";
        break;
      }
      iterations ++;
      x = aux;
      y = auy;
      d = FunctionUtils::derivative(f, x, y);

      if (verbose) {
        cout << "Iteration " << iterations << ": x = " << x << ", y = " << y << ", f'(x, y) = " << d << '\n';
      }

      if (iterations >= max_iters) {
        this->endReason = "Maximum number of iterations reached";
        break;
      }
      if (fabs(d) < error) {
        this->endReason = "Minimum error threshold reached";
        break;
      }
    }
    this->error = fabs(d);

    return {x, y};
  }
};

#endif //NUMERICAL_ANALYSIS_OPTIMIZER_HPP
