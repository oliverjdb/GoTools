#include <iostream>
#include <vector>

using namespace std;

// template<typename T> inline
// void operator += (T& t1, const T& t2) {
//   for (auto i1(t1.begin()), i2(t2.begin()); i1 != t1.end(); ++i1, ++i2)
//     *i1 += *i2;
// }
  
struct Verbose
{
  Verbose() : data(1) {cout << "Constructor called\n";}
  Verbose(const Verbose& rhs) : data(rhs.data) { cout << "Copy constructor called\n";}
  Verbose operator=(const Verbose& rhs) {
    cout << "Assignment operator called.\n";
    data = rhs.data;
  }
  
  int data;
};

Verbose operator += (Verbose& v1, const Verbose& v2) {
  v1.data += v2.data;
}

Verbose operator + (const Verbose& v1, const Verbose& v2)
{
  Verbose tmp;
  tmp.data = v1.data + v2.data;
  return tmp;
}

Verbose operator * (const Verbose& v, int i) {

  Verbose tmp;
  tmp.data = v.data * i;
  return tmp;
}
  

using namespace std;

int main() {

  Verbose v;
  Verbose v2(v);

  v += v2;

  cout << v.data <<endl;

  Verbose w = v + v2;

  cout << w.data <<endl;

  Verbose vv = (w + v2) * 3;

  cout << vv.data <<endl;
  
}

  
