#include <iostream>
#include <vector>

#define SILENCE_EXTERNAL_WARNINGS 1
#include "disable_warnings.h"
#include <Eigen/Dense>
#include "reenable_warnings.h"

using namespace std;
using namespace Eigen;

int main()
{
  //vector<double> v {1, 2, 3, 4, 5, 6, 7, 8, 9};
  Matrix<double, 3, 3> m(&(vector<double> { -1, -2, 3, 4, 5, 6, 7, 8, -9})[0]);
  std::cout << m << std::endl;
  std::cout << m.cwiseAbs().maxCoeff() << std::endl;

  return 0;
};
