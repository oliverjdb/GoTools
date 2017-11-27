#include <random>
#include <fstream>
#include "tritools.h"


using namespace std;
using namespace Go;
using namespace TriTools;

int main()
{
  // ---------------------- testing points_inside_triangle ----------------------
  Point p1 {0.0, 0.0};
  Point p2 {1.0, 0.0};
  Point p3 {0.0, 1.0};

  //generating random points
  size_t NUM_POINTS = 100000;
  random_device r;
  default_random_engine r_eng(r());
  uniform_real_distribution<double> uniform_dist(-1, 4);
  vector<double> point_coords(NUM_POINTS*2);
  generate(point_coords.begin(), point_coords.end(), [&]() {return uniform_dist(r_eng);});
  vector<Point> points(NUM_POINTS);

  // vector<Point> points(1);
  // points[0] = Point(0.1, 0.1);
  for (size_t i = 0; i != NUM_POINTS; ++i)
    points[i] = Point(point_coords[2*i], point_coords[2*i+1]);

  const auto points_inside = points_inside_triangle(p1, p2, p3, points, {0.3, 0.0, 0.0});

  // saving results
  {
    ofstream os1("all_points.dat");
    ofstream os2("inside_points.dat");
    for (size_t i = 0; i != NUM_POINTS; ++i) {
      os1 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
      if (points_inside[i])
	os2 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
    }
  }

  // ---------------------- testing boundary triangulation ----------------------
  
  vector<Point> loop_outer { {0.0, 0.0}, {3.0, 0.0}, {3.0,1.0},
  			     {2.0,1.0}, {2.0,2.0}, {3.0,2.0},
  			     {3.0,3.0}, {0.0,3.0}};

  vector<Point> loop_inner { {0.5, 1.0}, {0.5, 2.0}, {1.5, 2.0}, {1.5, 1.0}};

  // vector<Point> loop_outer { {0., 0.}, {3., 0.}, {3., 3.}, {0., 3.}};
  // vector<Point> loop_inner { {1., 1.}, {1., 2.}, {2., 2.}, {2., 1.}};
  
  vector<vector<Point>> loops;
  loops.push_back(loop_outer);
  loops.push_back(loop_inner);
  
  auto triangles = triangulate_boundary(loops);

  {
    ofstream os1("triang.dat");
    for (auto tri : triangles) {
      os1 << tri.first[0] << " " << tri.first[1] << " " << tri.first[2] << '\n';
    }
  };
    
  // ------------------------ Testing points inside loops ------------------------
  
  auto ind = points_inside_loops(loops, points, 0.2);
  {
    ofstream os1("all_points2.dat");
    ofstream os2("inside_points2.dat");
    for (size_t i = 0; i != NUM_POINTS; ++i) {
      os1 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
      if (ind[i])
	os2 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
    }
  }
    
  return 0;
}
