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

  double s1 = expm1(1.0), s2 = M_PI_4, s3 = sqrt(M_PI) / 2 * erf(high);

  cout << expm1(1.0) << endl;
  cout << o.integrate(fe, low, high, quadratures, Optimizer::RECTANGLE) << endl;
  cout << o.integrate(fe, low, high, quadratures, Optimizer::TRAPEZOID) << endl;
  cout << o.integrate(fe, low, high, quadratures, Optimizer::SIMPSON) << endl;
  try {
    cout << o.adaptiveIntegration(fe, low, high, Optimizer::RECTANGLE, s1) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(fe, low, high, Optimizer::TRAPEZOID, s1) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(fe, low, high, Optimizer::SIMPSON, s1) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }

  cout << s2 << endl;
  cout << o.integrate(ff, low, high, quadratures, Optimizer::RECTANGLE) << endl;
  cout << o.integrate(ff, low, high, quadratures, Optimizer::TRAPEZOID) << endl;
  cout << o.integrate(ff, low, high, quadratures, Optimizer::SIMPSON) << endl;
  try {
    cout << o.adaptiveIntegration(ff, low, high, Optimizer::RECTANGLE, s2) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(ff, low, high, Optimizer::TRAPEZOID, s2) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(ff, low, high, Optimizer::SIMPSON, s2) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }

  cout << s3 << endl;
  cout << o.integrate(fg, low, high, quadratures, Optimizer::RECTANGLE) << endl;
  cout << o.integrate(fg, low, high, quadratures, Optimizer::TRAPEZOID) << endl;
  cout << o.integrate(fg, low, high, quadratures, Optimizer::SIMPSON) << endl;
  try {
    cout << o.adaptiveIntegration(fg, low, high, Optimizer::RECTANGLE, s3) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(fg, low, high, Optimizer::TRAPEZOID, s3) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(fg, low, high, Optimizer::SIMPSON, s3) << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }

  return 0;
}
