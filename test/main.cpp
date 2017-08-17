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

//! \param x
//! \return e^x
double fe(double x) { return exp(x); }

//! \param x
//! \return sqrt(1 - pow(x, 2))
double ff(double x) { return sqrt(1 - pow(x, 2)); }

//! \param x
//! \return exp(-x^2)
double fg(double x) { return exp(-(x*x)); }

void testRoots(double x, double error, double result, int iters,
               tuple<double, double> &d2_ret, int learnRateFraction,
               Optimizer &o) {
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.findRoot(fa, x, error, iters, learnRate);
    cout << "Root,1," << learnRate << "," << result << "," << o.getIterations()
         << "," << o.getError() << "," << o.getEndReason() << "\n";
  }
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.findRoot(fb, x, error, iters, learnRate);
    cout << "Root,2," << learnRate << "," << result << "," << o.getIterations()
         << "," << o.getError() << "," << o.getEndReason() << "\n";
  }
}

void testMinimization(double x, double y, double error, double result,
                      int iters, tuple<double, double> &d2_ret,
                      int learnRateFraction, Optimizer &o) {
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.minimize(fa, x, error, iters, learnRate);
    cout << "Minimum,1," << learnRate << "," << result << ","
         << o.getIterations() << "," << o.getError() << "," << o.getEndReason()
         << "\n";
  }
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    result = o.minimize(fb, x, error, iters, learnRate);
    cout << "Minimum,2," << learnRate << "," << result << ","
         << o.getIterations() << "," << o.getError() << "," << o.getEndReason()
         << "\n";
  }
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    d2_ret = o.minimize(fc, x, y, error, iters, learnRate);
    cout << "Minimum,3," << learnRate << "," << get<0>(d2_ret) << ", "
         << get<1>(d2_ret) << "," << o.getIterations() << "," << o.getError()
         << "," << o.getEndReason() << "\n";
  }
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    d2_ret = o.minimize(fd, x, y, error, iters, learnRate);
    cout << "Minimum,4," << learnRate << "," << get<0>(d2_ret) << ", "
         << get<1>(d2_ret) << "," << o.getIterations() << "," << o.getError()
         << "," << o.getEndReason() << "\n";
  }
}

int main() {
  cout.precision(12);
  double x = 2, y = 2, error = 1e-50, result, low = 0, high = 1;
  int iters = 1000000;
  tuple<double, double> d2_ret;
  int learnRateFraction = 100;
  int quadratures = 200;
  Optimizer o;

  //  testRoots(x, error, result, iters, d2_ret, learnRateFraction, o);
  //  testMinimization(x, y, error, result, iters, d2_ret, learnRateFraction, o);

  }

  return 0;
}
