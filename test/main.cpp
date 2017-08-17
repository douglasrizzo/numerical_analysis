#include "FunctionUtils.hpp"
#include "Optimizer.hpp"

using namespace std;

//! \param x
//! \return x^2
double fa(double x) { return x * x; }

//! \param x
//! \return x ^ 3 - 2x ^ 2 + 2
double fb(double x) { return pow(x, 3) - 2 * pow(x, 2) + 2; }

//! \param x
//! \param y
//! \return (1 - x) ^ 2 + (1 - y) ^ 2
double fc(double x, double y) { return pow((1 - x), 2) + pow((1 - y), 2); }

//! Calculates Rosenbrock's function
//! \param x
//! \param y
//! \return (1 - y) ^ 2 + 100(x - y ^ 2) ^ 2
double fd(double x, double y) {
  return pow((1 - y), 2) + 100 * pow((x - pow(y, 2)), 2);
}

int main() {
  cout.precision(12);
  double x = 2, y = 2, error = 1e-8, result;
  int iters = 1000;
  tuple<double, double> d2_ret;
  int learnRateFraction = 100;
  Optimizer o;

  cout << "Finding roots of x^2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.findRoot(fa, x, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << result
         << " (0), " << o.getIterations() << " iterations, error = " << o.getError() << ", end reason: "
         << o.getEndReason() << "\n";
  }
  cout << "\nFinding roots of x^3 - 2x^2 + 2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.findRoot(fb, x, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << result
         << " (-0.83929), " << o.getIterations() << " iterations, error = " << o.getError() << ", end reason: "
         << o.getEndReason() << "\n";
  }
  cout << "\nFinding minimum of x^2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.minimize(fa, x, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << result
         << " (0), " << o.getIterations() << " iterations, error = " << o.getError() << ", end reason: "
         << o.getEndReason() << "\n";
  }
  cout << "\nFinding minimum of x^3 - 2x^2 + 2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.findRoot(fa, x, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << result
         << " (" << 4.0 / 3.0 << "), " << o.getIterations() << " iterations, error = " << o.getError()
         << ", end reason: " << o.getEndReason() << "\n";
  }
  cout << "\nFinding minimum of (1 - x)^2 + (1 - y)^2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    d2_ret = o.minimize(fc, x, y, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << '(' << get<0>(d2_ret) << ", " << get<1>(d2_ret) << ") (1, 1), "
         << o.getIterations() << " iterations, error = " << o.getError() << ", end reason: " << o.getEndReason()
         << "\n";
  }
  cout << "\nFinding minimum of (1 - y)^2 + 100(x - y^2)^2...\n";
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    d2_ret = o.minimize(fd, x, y, error, iters, learnRate);
    cout << "Beta = " << learnRate << ": x = " << '(' << get<0>(d2_ret) << ", " << get<1>(d2_ret)
         << ") (1, 1), " << o.getIterations() << " iterations, error = " << o.getError() << ", end reason: "
         << o.getEndReason() << "\n";
  }

  return 0;
}
