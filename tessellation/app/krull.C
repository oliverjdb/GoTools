#include <fstream>
#include <iostream>
#include <string>
#include "tesselate_curve.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace std;
using namespace Go;

namespace {
  void read_away_header(ifstream& is) {
    for (int i=0, tmp=0;i != 4; is >> tmp, ++i); // read away header  
  }
  
}; // anonymous namespace 

int main(int varnum, char* vararg[]) {

  // ------------------------------ loading object ------------------------------
  const string filename = "krull2.g2";

  ifstream is(filename.c_str());

  vector<SplineCurve> curves(4);
  
  for (int i = 0; i != 4; ++i) {
    read_away_header(is);
    curves[i].read(is);
  }

  vector<SplineCurve> flat_curves = curves;
  for (int i = 0; i != 4; ++i) {
    SplineCurve& s = flat_curves[i];
    for (int j = 0; j != s.numCoefs(); ++j)
      *(s.coefs_begin() + 3*j + 2) = 0;
  }

  // ofstream os("krull2_flat_mollified.g2");
  // for (int i = 0; i != 4; ++i) {
  //   flat_curves[i].writeStandardHeader(os);
  //   flat_curves[i].write(os);
  // }
  // os.close();

  ifstream is2("krull2_flat_mollified.g2");
  for (int i = 0; i != 4; ++i) {
    read_away_header(is2);
    flat_curves[i].read(is2);
  }
  
  //tesselate_curve(flat_curves[1], 20);
  tesselate_curve(flat_curves[3], (unsigned int)10);



  // for (auto c : curve_loop) {
  //   auto krull = c->geometryCurve();
  //   krull->writeStandardHeader(os);
  //   krull->write(os);
  //   //tesselate_curve(*(c->geometryCurve()), 15);
  // }
  // os.close();
  

  
  return 0;
}
  
  
