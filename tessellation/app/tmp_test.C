#include <array>
#include <iostream>

using namespace std;

void operator+=(array<double, 4>& a1, const array<double,4>& a2) 
{
  auto it2 = a2.cbegin();
  for(auto it = a1.begin(); it != a1.end(); ++it, ++it2)
    *it += *it2;
}

int main() {

  const array<double, 4> a1 {1, 2, 3, 4};
  const array<double, 4> a2 {5, 4, 0, 2};
  
  auto result = a1;
  result += a2;

  for (auto v : result)
    cout << v << " ";


  return 0;
};
