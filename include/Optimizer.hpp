/**
 * @author Douglas De Rizzo Meneghetti (douglasrizzom@gmail.com)
 * @brief  Numerical optimizer specialized in finding roots, minima and integrals of functions
 * @date   2017-8-11
 */

#ifndef NUMERICAL_ANALYSIS_OPTIMIZER_HPP
#define NUMERICAL_ANALYSIS_OPTIMIZER_HPP

#include "FunctionUtils.hpp"
#include <cmath>
#include <functional>
#include <iostream>

//! Numerical optimizer specialized in finding roots, minima and integrals of functions
class Optimizer {
 public:
  enum IntegrationMethod { RECTANGLE, TRAPEZOID, SIMPSON };

 private:
  int iterations = 0;
  double error = 0;
  string endReason = "You didn't run any optimization yet!";


  //! Adaptive quadrature recursive method
  //! \param f function to integrate
  //! \param a the lower bound of the integration interval
  //! \param b the upper bound of the integration interval
  //! \param method the Newton-Cotes function to use in the approximation
  //! \param error
  //! \return Numerical approximation of the integral of f
  double innerAdaptiveIntegration(const function<double(double)> &f, double a,
                                  double b, IntegrationMethod method,
                                  double error = 1e-12) {
    iterations += 2;
    // calculates the middle point between a and b
    double meio = (b + a) / 2;
    // uses the integrate method to calculate the value of a single quadrature
    // vs. the sum of two sub-quadratures
    double i1 = integrate(f, a, b, 1, method),
        i2 = integrate(f, a, meio, 1, method) +
        integrate(f, meio, b, 1, method);

    // if there is error, run adaptive integration in the two sub-divisions of
    // the current partition
    if (fabs(i1 - i2) > error) {
      return innerAdaptiveIntegration(f, a, meio, method, error) +
          innerAdaptiveIntegration(f, meio, b, method, error);
    }

    // otherwise, return the most precise value of the two already calculated
    return i2;
  }

 public:
  //! \return number of iterations until convergence
  int getIterations() const { return iterations; }

  //! \return error of approximation
  double getError() const { return error; }

  //! \return reason for ending the process
  const string &getEndReason() const { return endReason; }

  //! Numerically searches for the root of a function via Newton-Raphson method
  //! \param f a function
  //! \param x the initial point to start the search
  //! \param error minimum tolerance for the search to end
  //! \param max_iters maximum number of iterations
  //! \param learnRate the learning rate of the search
  //! \param verbose whether to print a short summary of the search at every iteration
  //! \return the point at which the function intercepts the x-axis
  double findRoot(const std::function<double(double)> &f, double x,
                  double error = 1e-8, int max_iters = 1000,
                  double learnRate = 1, bool verbose = false) throw(runtime_error) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double f_val;

    while (true) {
      f_val = f(x);
      double aux = x + learnRate * - f_val / FunctionUtils::derivative(f, x);
      if (aux == x) {
        this->endReason = "No change in x from previous iteration";
        break;
      }

      iterations ++;
      x = aux;
      if (verbose) {
        cout << "Iteration " << iterations << ": x = " << x
             << ", f(x) = " << f_val << '\n';
      }

      if (iterations >= max_iters) {
        string message = "Maximum number of iterations reached";
        this->endReason = message;
        throw runtime_error(message);
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
                  double error = 1e-8, int max_iters = 1000,
                  double learnRate = 1, bool verbose = false) throw(runtime_error) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double d;

    while (true) {
      d = FunctionUtils::derivative(f, x);
      double aux = x - learnRate * d;
      if (aux == x) {
        this->endReason = "No change in x from previous iteration";
        break;
      }
      iterations ++;
      x = aux;

      if (verbose) {
        cout << "Iteration " << iterations << ": x = " << x << ", f'(x) = " << d
             << '\n';
      }

      if (iterations >= max_iters) {
        string message = "Maximum number of iterations reached";
        this->endReason = message;
        throw runtime_error(message);
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
  //! \param verbose whether to print a short summary of the search at every
  //! iteration
  //! \return a tuple containing the {x, y} points at which the function is
  //! minimal
  std::tuple<double, double>
  minimize(const std::function<double(double, double)> &f, double x, double y,
           double error = 1e-8, int max_iters = 1000, double learnRate = 1,
           bool verbose = false) throw(runtime_error) {
    endReason = "You didn't run any optimization yet!";
    iterations = 0;
    double dfdx, dfdy;

    while (true) {
      dfdx = FunctionUtils::partialDerivative(f, x, y, 0);
      dfdy = FunctionUtils::partialDerivative(f, x, y, 1);

      double aux = x - learnRate * dfdx;
      double auy = y - learnRate * dfdy;

      if (aux == x and y == auy) {
        this->endReason = "No change in x and y from previous iteration";
        break;
      }

      iterations ++;

      x = aux;
      y = auy;

      if (verbose and iterations % 1000000 == 0) {
        cout << "x = " << x << ", y = " << y << ", f'(x, y) = (" << dfdx << ", "
             << dfdy << ")\tIteration " << iterations << endl;
      }

      if (iterations >= max_iters) {
        string message = "Maximum number of iterations reached";
        this->endReason = message;
        throw runtime_error(message);
      }
      if (fabs(dfdx) + fabs(dfdy) < error) {
        this->endReason = "Minimum error threshold reached";
        break;
      }
    }
    this->error = (fabs(dfdx) + fabs(dfdy) / 2);

    return {x, y};
  }

  //! Numerically approximates the integral of a function
  //! \param f the function to integrate
  //! \param low the lower bound of the integration interval
  //! \param high the upper bound of the integration interval
  //! \param points the number of quadrature points to use in the approximation
  //! \param method the Newton-Cotes function to use in the approximation
  //! \return Numerical approximation of the integral of f
  double integrate(const std::function<double(double)> &f, double low,
                   double high, long int points = 40,
                   IntegrationMethod method = SIMPSON) throw(runtime_error) {

    if (low == high) {
      throw runtime_error("Lower bound of integration = Higher bound");
    }
    if (low > high) {
      double temp = low;
      low = high;
      high = temp;
    }

    std::function<double(std::function<double(double)>, double, double)> approx;

    if (method == SIMPSON)
      approx = FunctionUtils::simpsonRule;
    else if (method == RECTANGLE)
      approx = FunctionUtils::rectangleRule;
    else if (method == TRAPEZOID)
      approx = FunctionUtils::trapezoidRule;
    else
      throw runtime_error("Unsupported integration method");

    double sum = 0;
    double step = (high - low) / points;

    if (step < 1e-8)
      throw runtime_error("Step size of " + to_string(step) + " is too small to be precise");

    for (int i = 0; i < points; i ++)
      sum += approx(f, low + (i * step), low + ((i + 1) * step));

    return sum;
  }

  //! Adaptive quadrature method
  //! \param f function to integrate
  //! \param a the lower bound of the integration interval
  //! \param b the upper bound of the integration interval
  //! \param method the Newton-Cotes function to use in the approximation
  //! \param error
  //! \return Numerical approximation of the integral of f
  double adaptiveIntegration(const function<double(double)> &f, double a,
                             double b, IntegrationMethod method,
                             double error = 1e-12) {
    iterations = 0;
    return innerAdaptiveIntegration(f, a, b, method, error);
  }
};

#endif // NUMERICAL_ANALYSIS_OPTIMIZER_HPP
