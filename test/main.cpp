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
double fg(double x) { return exp(- (x * x)); }

void testSingleRoot(const function<double(double)> &f,
                    double x,
                    double error,
                    int iters,
                    int learnRateFraction) {
  Optimizer o;
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;
    double result = o.findRoot(f, x, error, iters, learnRate);
    cout << "Root,1," << learnRate << "," << result << "," << o.getIterations()
         << "," << o.getError() << "," << o.getEndReason() << "\n";
  }
}

void testRoots(double x, double error, int iters,
               int learnRateFraction,
               Optimizer &o) {
  testSingleRoot(fa, x, error, iters, learnRateFraction);
  testSingleRoot(fb, x, error, iters, learnRateFraction);
}

void testSingleVariableMinimization(const function<double(double)> &f,
                                    double x,
                                    double error,
                                    int iters,
                                    int learnRateFraction) {
  Optimizer o;
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;

    try {
      double result = o.minimize(f, x, error, iters, learnRate);
      cout << "Minimum,1," << learnRate << "," << result << ","
           << o.getIterations() << "," << o.getError() << "," << o.getEndReason()
           << "\n";
    } catch (const runtime_error &exp) {
      cout << exp.what() << endl;
    }
  }
}

void testDoubleVariableMinimization(const function<double(double, double)> &f, double x,
                                    double y,
                                    double error,
                                    int iters,
                                    int learnRateFraction) {
  Optimizer o;
  tuple<double, double> d2_ret;
  for (int i = 1; i <= learnRateFraction; i ++) {
    double learnRate = (double) i / learnRateFraction;

    try {
      d2_ret = o.minimize(f, x, y, error, iters, learnRate);
      cout << "Minimum,3," << learnRate << "," << get<0>(d2_ret) << ", "
           << get<1>(d2_ret) << "," << o.getIterations() << "," << o.getError()
           << "," << o.getEndReason() << "\n";
    } catch (const runtime_error &exp) {
      cout << exp.what() << endl;
    }
  }
}

void testMinimization(double x, double y, double error,
                      int iters,
                      int learnRateFraction) {
  testSingleVariableMinimization(fa, x, error, iters, learnRateFraction);
  testSingleVariableMinimization(fb, x, error, iters, learnRateFraction);
  testDoubleVariableMinimization(fc, x, y, error, iters, learnRateFraction);
  testDoubleVariableMinimization(fd, x, y, error, iters, learnRateFraction);
}

void testSingleIntegral(const function<double(double)> &f, double low, double high, int quadratures, double trueValue) {
  Optimizer o;
  cout << trueValue << "\t\"true\" value" << endl;
  cout << o.integrate(f, low, high, quadratures, Optimizer::RECTANGLE) << "\trectangle rule" << endl;
  cout << o.integrate(f, low, high, quadratures, Optimizer::TRAPEZOID) << "\ttrapezoid rule" << endl;
  cout << o.integrate(f, low, high, quadratures, Optimizer::SIMPSON) << "\tsimpson rule" << endl;
  try {
    cout << o.adaptiveIntegration(f, low, high, Optimizer::RECTANGLE, trueValue) << "\tadaptive rectangle rule" << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(f, low, high, Optimizer::TRAPEZOID, trueValue) << "\tadaptive trapezoid rule" << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
  try {
    cout << o.adaptiveIntegration(f, low, high, Optimizer::SIMPSON, trueValue) << "\tadaptive simpson rule" << endl;
  } catch (const runtime_error &x) {
    cout << x.what() << endl;
  }
}

void testIntegrals(double low, double high, int quadratures) {
  double s1 = expm1(1.0), s2 = M_PI_4, s3 = sqrt(M_PI) / 2 * erf(high);

  cout << "Integrating e^x..." << endl;
  testSingleIntegral(fe, low, high, quadratures, s1);
  cout << "Integrating sqrt(1 - pow(x, 2))..." << endl;
  testSingleIntegral(ff, low, high, quadratures, s2);
  cout << "Integrating exp(-(x^2))..." << endl;
  testSingleIntegral(fg, low, high, quadratures, s3);
}

int main() {
  cout.precision(12);
  double x = 2, y = 2, error = 1e-50, low = 0, high = 1;
  int iters = 1000000;
  int learnRateFraction = 100;
  int quadratures = 200;
  Optimizer o;

  testRoots(x, error, iters, learnRateFraction, o);
  testMinimization(x, y, error, iters, learnRateFraction);
  testIntegrals(low, high, quadratures);

  return 0;
}
