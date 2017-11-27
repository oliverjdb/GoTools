#include <fstream>
#include <iostream>
#include <string>
#include <iterator>
#include "tesselate_curve.h"
#include "tesselate_surface.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"

#include "tesselate_debug.h"
#include "tritools.h"
#include "make_spacing_funs.h"

using namespace std;
using namespace Go;
using namespace TriTools;

namespace {
  void read_away_header(ifstream& is) {
    for (int i=0, tmp=0;i != 4; is >> tmp, ++i); // read away header  
  }

  // const vector<Point> construct_gridpoints(double minx, double maxx, size_t samples_x,
  // 					   double miny, double maxy, size_t samples_y)
  // {
  //   vector<Point> result;
  //   for (size_t iy = 0; iy != samples_y; ++iy)
  //     for (size_t ix = 0; ix != samples_x; ++ix)
  // 	result.push_back(Point {minx + (maxx - minx)/double(samples_x - 1) * double(ix),
  // 	                        miny + (maxy - miny)/double(samples_y - 1) * double(iy)});
  //   return result;
  // }
  
}; // anonymous namespace 

int main(int varnum, char* vararg[]) {

  // ------------------------------ loading object ------------------------------
  const string filename = (varnum > 1) ? vararg[1] : "curve_and_surface.g2";

  ifstream is(filename.c_str());

  read_away_header(is);
  shared_ptr<ParamSurface> ss {new SplineSurface};
  ss->read(is);

  // read four curves

  GoTools::init();
  vector<shared_ptr<CurveOnSurface>> curve_loop;
  for (size_t i = 0; i != 4; ++i) {
    read_away_header(is);
    shared_ptr<CurveOnSurface> c {new CurveOnSurface};
    c->read(is);
    c->setUnderlyingSurface(ss);
    curve_loop.emplace_back(c);
  }

  // ----------------------------- tesselate curves -----------------------------

  vector<Point> param_points; // points describing the closed curve polygon in
			      // parameter space
  ofstream os("parampoints.dat");
  for (auto c : curve_loop) {
    //auto krull = c->geometryCurve();
    //vector<double> pvec = tesselate_curve(*(c->geometryCurve()), (unsigned int)20);
    vector<double> pvec = tesselate_curve(*(c->geometryCurve()), 0.1);
    
    // computing points
    Point pt;
    c->parameterCurve()->point(pt, c->startparam());
    param_points.push_back(pt);
    for (auto p = pvec.begin(); p != pvec.end(); ++p) {
      c->parameterCurve()->point(pt, *p);
      param_points.push_back(pt);
    }
    //store_points_and_curve(*krull, &pvec[0], pvec.size(), "dill.g2");
  }
  for (auto p : param_points)
    os << p[0] << " " << p[1] << endl;
  
  os.close();

  // ----- choose interior sample locations and clip against bounding curve -----

  vector<vector<Point>> loops;
  loops.push_back(param_points);
  SurfaceTriangulation tri = tesselate_surface(*ss, 25, 25, 0.1, loops);
  
  // const size_t SAMPLES = 50;
  // const vector<Point> gridpoints = construct_gridpoints(0, 1, SAMPLES, 0, 1, SAMPLES);
  // const CurveBoundedDomain domain = construct_curve_bounded_domain(param_points);
  // vector<Point> kept_points;
  // for (p : gridpoints)
  //   if (domain.isInDomain(p, 1e-6))
  //     kept_points.push_back(p);

  // ofstream os2("intpoints.dat");
  // for (auto p : kept_points)
  //   os2 >> p[0] << " " << p[1] << endl;

  ofstream os3("triangles.dat");
  for (auto t : tri.triangles) {
    os3 << tri.uv[t[0]][0] << " " << tri.uv[t[0]][1] << " ";
    os3 << tri.uv[t[1]][0] << " " << tri.uv[t[1]][1] << " ";
    os3 << tri.uv[t[2]][0] << " " << tri.uv[t[2]][1] << "\n";
  }


  vector<double> par(tri.num_interior_nodes * 2, 0l);
  for (size_t i = 0; i != tri.num_interior_nodes; ++i) {
    par[2*i]   = tri.uv[i][0];
    par[2*i+1] = tri.uv[i][1];
  }
  
  // vector<Point> krullpoints {Point {1.0, 1.0}, {0.0,0.0}, {2.0,0.0}, {2.0, 2.0}, {0.0, 2.0}};
  // vector<array<int,3>> krulltris {{0,1,2}, {0,2,3}, {0,3,4}, {0,4,1}};
  // SurfaceTriangulation krull {krullpoints, krulltris, (size_t)1};
  // vector<double> krullpar {1.1, 1.1};

  // auto tsfun = make_triang_spacing_fun(krull);
  // auto funval = get<0>(tsfun)(&krullpar[0], 0);
  // auto dfunval = get<1>(tsfun)(&krullpar[0], 0);
  auto tsfun = make_triang_spacing_fun(tri);
  auto funval = get<0>(tsfun)(&par[0], 0);
  auto dfunval = get<1>(tsfun)(&par[0], 0);
  ofstream dfun1("gradient.dat");
  ofstream dfun2("jacobian.dat");
  copy(dfunval.value.begin(), dfunval.value.end(), ostream_iterator<double>(dfun1, " "));
  copy(dfunval.jacobian.begin(), dfunval.jacobian.end(), ostream_iterator<double>(dfun2, " "));
  return 0;
}
  
  
