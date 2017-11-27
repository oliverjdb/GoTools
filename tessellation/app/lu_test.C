#include <iterator>
#include <algorithm>
#include <iostream>
#include <vector>
#include "lu.h"


using namespace std;
using namespace TesselateUtils;

int main() {

  vector<double> coefs {1, 2, 3, 4, 5, 6 ,7, 8, 9};
  vector<int> perm(3, 0);
  bool parity;

  lu(3, &coefs[0], &perm[0], parity);

  for (unsigned int i = 0; i != 3; ++i, cout << '\n')
    for (unsigned int j = 0; j != 3; ++j)
      cout << coefs[i + j*3] << " ";

  cout << endl;
  copy(perm.begin(), perm.end(), ostream_iterator<int>(cout, " "));

  return 0;
}

  
