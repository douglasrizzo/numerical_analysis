#include <iostream>
#include <cmath>
#include <FunctionUtils.hpp>

using namespace std;

double fa(double x) {
  return pow(x, 2);
}

double fb(double x) {
  return pow(x, 3) - 2 * pow(x, 2) + 2;
}

double fc(double x, double y) {
  return pow((1 - x), 2) + pow((1 - y), 2);
}

double fd(double x, double y) {
  return pow((1 - y), 2) + 100 * pow((x - pow(y, 2)), 2);
}

int main() {

  double x = 7, y = 7, delta = 0.03;
  std::tuple<double, double> d2_ret;
  for (int i = 1; i <= 10; i ++) {
    double learnRate = (double) i / 10;
    cout << learnRate << endl;
    cout << FunctionUtils::findRoot(fa, x, delta, 1e-5, learnRate) << " (0)\n";
    cout << FunctionUtils::findRoot(fb, x, delta, 1e-5, learnRate) << " (-0.839287)\n";
    cout << FunctionUtils::minimize(fa, x, delta, 1e-5, learnRate) << " (0)\n";
//    cout << FunctionUtils::minimize(fb, x, delta, 1e-5, learnRate) << " (1.5)\n";
    d2_ret = FunctionUtils::minimize(fc, x, y, delta, 1e-5, learnRate);
    cout << '(' << std::get<0>(d2_ret) << ", " << std::get<1>(d2_ret) << ") (1, 1)\n";
    d2_ret = FunctionUtils::minimize(fd, x, y, delta, 1e-5, learnRate);
    cout << '(' << std::get<0>(d2_ret) << ", " << std::get<1>(d2_ret) << ") (1, 1)\n\n";
  }
  return 0;
}
