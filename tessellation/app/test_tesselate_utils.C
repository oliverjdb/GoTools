#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include "tesselate_utils.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/CPUclock.h"
#include "interpoint_distances.h"
#include "polyhedral_energies.h"
#include "triangulate_domain.h"
#include "SimplePolyhedronTesselation.h"
#include "GoParametricTesselableVolume.h"
#include "tesselate_polyhedron.h"
#include "fit_points_to_plane.h"
#include "basic_intersections.h"
#include "clip_grid.h"
#include "debug_tesselate.h"

#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
//#include "GoTools/compositemodel/CompositeModelFileHandler.h"
#include "GoTools/compositemodel/SurfaceModel.h"

using namespace Go;
using namespace std;
using namespace TesselateUtils;


namespace Test { 
  void test_poly_p();
  void test_bounding_box_2D();
  void test_generate_grid();
  void test_inpolygon();
  void test_interpoint_distances();
  void test_energy();
  void test_circumscribe_triangle();
  void test_segment_intersection();
  void test_triangulation();
  void test_volume_tesselation();
  void test_fit_to_plane();
  void test_linear_system();
  void test_polyhedron3D_tesselation();
  void test_circumsphere();
  void test_triangle_intersection();
  void test_solve_linear_system();
  void test_clip_grid_2D();
  void test_clip_grid_3D();
  void test_parametric_volume_tesselation();
};

// ============================================================================
int main() {

  // cout << "Testing polygon area: " << endl;
  // Test::test_poly_area();

  // cout << "Testing bounding box: " << endl;
  // Test::test_bounding_box_2D();

  // cout << "Testing generate_grid: " << endl;
  // Test::test_generate_grid();

  // cout << "Testing inpolygon: " << endl;
  // Test::test_inpolygon();
  
  // cout << "Testing interpoint distance computation: " << endl;
  // Test::test_interpoint_distances();

  // cout << "Testing energy computation: " << endl;
  // Test::test_energy();

  // cout << "testing circumscribe triangle: " << endl;
  // Test::test_circumscribe_triangle();

  // cout << "testing segment intersection: " << endl;
  // Test::test_segment_intersection();
  
  // cout << "testing triangulation generation: " << endl;
  // Test::test_triangulation();

  // cout << "testing plane fitting: " << endl;
  // Test::test_fit_to_plane();

  // cout << "testing volume tesselation: " << endl;
  // Test::test_volume_tesselation();

  cout << "testing parametric volume tesselation: " << endl;
  Test::test_parametric_volume_tesselation();
  
  // cout << "testing linear system solving: " << endl;
  // Test::test_linear_system();

  // cout << "Testing tesselation of 3D polyhedron: " << endl;
  // Test::test_polyhedron3D_tesselation();

  // cout << "Testing circumsphere: " << endl;
  // Test::test_circumsphere();
  
  // Test::test_triangle_intersection();

  //Test::test_solve_linear_system();

  //Test::test_clip_grid_2D();

  //Test::test_clip_grid_3D();
  
  return 0;
};

// ============================================================================

namespace Test {
  
// ----------------------------------------------------------------------------
void test_poly_area()
// ----------------------------------------------------------------------------
{
  vector<Go::Point> poly { {0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}};

  cout << "Area is: " << polygon_area(&poly[0], (int)poly.size()) << endl;
}

// ----------------------------------------------------------------------------
void test_bounding_box_2D()
// ----------------------------------------------------------------------------
{
  vector<Go::Point> poly { {0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}};
  
  const auto box = bounding_box_2D(&poly[0], (int)poly.size());
  cout << "Bounding box is: \n";
  cout << "xmin: " << box[0] << ", xmax: " << box[1] << '\n';
  cout << "ymin: " << box[2] << ", ymax: " << box[3] << endl;
}

// ----------------------------------------------------------------------------  
void test_generate_grid()
// ----------------------------------------------------------------------------  
{
  const Go::Point p1 {0.0, 0.0};
  const Go::Point p2 {1.0, 2.0};
  const int nx = 4;
  const int ny = 6;

  const auto result = generate_grid_2D(p1, p2, nx, ny);

  cout << "Gridpoints are: \n";
  copy(result.begin(), result.end(), std::ostream_iterator<Go::Point>(cout, "\n"));
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_inpolygon()
// ----------------------------------------------------------------------------  
{
  const vector<Go::Point> poly = { {0.0, 0.0}, {1.0, 0.0}, 
				   {1.0, 2.0}, {0.5, 1.0}, {0.0, 1.0}};
  const auto bbox = bounding_box_2D(&poly[0], (int)poly.size());
  const int nx = 30;
  const int ny = 30;
  const double tol = 1e-6;
  const auto candidates = generate_grid_2D(Go::Point {bbox[0], bbox[2]}, 
					   Go::Point {bbox[1], bbox[3]}, nx, ny);
  const auto result = inpolygon(&candidates[0], (unsigned int)candidates.size(),
				&poly[0], (unsigned int)poly.size(), tol);

  cout << "Polygon is: \n";
  copy(poly.begin(), poly.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

  cout << "Candidate points were: \n";
  copy(candidates.begin(), candidates.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

  cout << "Interior points are: \n" << endl;
  copy(result.begin(), result.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

}

// ----------------------------------------------------------------------------  
void test_interpoint_distances()
// ----------------------------------------------------------------------------
{
  // Define a random set of test points
  const int N = 20000;
  const double R = 0.1;
  vector<Point2D> points(N);
  generate(points.begin(), points.end(), 
	   [](){return Point2D {random_uniform(0, 1), random_uniform(0, 1)};});

  const auto result = interpoint_distances(&points[0], N, R, false);

  // ofstream os("distances.mat");
  // for (auto e : result)
  //   os << e.p1_ix << " " << e.p2_ix << " " << e.dist << '\n';
  // os.close();

  // ofstream os2("points.mat");
  // for (auto p : points)
  //   os2 << p[0] << " " << p[1] << '\n';
  // os2.close();

  cout << "Number of relations found: " << result.size() << endl;
    
}

// ----------------------------------------------------------------------------
void test_energy()
// ----------------------------------------------------------------------------
{
  const vector<Point2D> bpoints { {0,0}, {1,0}, {1,1}, {0,1} };

  const vector<Point2D> points { {0.25, 0.5}, {0.75, 0.5}};
  const double R = 1;

  const auto result = polygon_energy(&bpoints[0], (uint)bpoints.size(),
				     &points[0],  (uint)points.size(),
				     R);

  cout << "Energy is: " << result.val << '\n';

  cout << "Energy derivatives are: \n";
  for (auto d : result.der) {
    cout << d[0] << " " << d[1] << '\n';
  }
    
}
// ----------------------------------------------------------------------------
void test_circumscribe_triangle()
// ----------------------------------------------------------------------------
{
  const Point2D p1 {0, 0};
  const Point2D p2 {1, 0};
  const Point2D p3 {1, 1};

  Point2D center;
  double radius2;
  circumscribe_triangle(p1, p2, p3, center, radius2);

  cout << "Points are: \n";
  cout << "(" << p1[0] << ", " << p1[1] << "), ";
  cout << "(" << p2[0] << ", " << p2[1] << "), ";
  cout << "(" << p3[0] << ", " << p3[1] << ")\n";
  cout << "Center is: \n";
  cout << "(" << center[0] << ", " << center[1] << ")\n";
  cout << "Radius is: " << sqrt(radius2) << endl << endl;

  cout << "Computed distances to each point are: \n";
  cout << "For p1 -- " << dist(p1, center) << endl;
  cout << "For p2 -- " << dist(p2, center) << endl;
  cout << "For p3 -- " << dist(p3, center) << endl;
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_segment_intersection()
// ----------------------------------------------------------------------------
{
  // obvious case
  const Point2D p1a {1, 1}, p1b {4, 3};
  const Point2D p2a {3, 1}, p2b {2, 3};

  // endpoints-meet case
  const Point2D q1a {1, 1}, q1b {1, 4};
  const Point2D q2a {1, 4}, q2b {2, 3};

  // endpoint-touch-segment case
  const Point2D r1a {1, 1}, r1b {3, 1};
  const Point2D r2a {2, 1}, r2b {2, 2};

  const double tol = 1e-6;
  cout << "Obvious case, positive tol: " << segments_intersect_2D(p1a, p1b, p2a, p2b,  tol) << endl;
  cout << "Obvious case, negative tol: " << segments_intersect_2D(p1a, p1b, p2a, p2b, -tol) << endl;
  cout << endl;
  cout << "Endpoint-meet, positive tol: " << segments_intersect_2D(q1a, q1b, q2a, q2b,  tol) << endl;
  cout << "Endpoint-meet, negative tol: " << segments_intersect_2D(q1a, q1b, q2a, q2b, -tol) << endl;
  cout << endl;
  cout << "Point-touch-line, pos. tol: " << segments_intersect_2D(r1a, r1b, r2a, r2b,  tol) << endl;
  cout << "Point-touch-line, pos. tol: " << segments_intersect_2D(r1a, r1b, r2a, r2b, -tol) << endl;
  cout << endl;
  
  
}

// ----------------------------------------------------------------------------
void test_triangulation()
// ----------------------------------------------------------------------------
{
  vector<Point2D> bpoints { {0, 0}, {1, 0}, {1, 1}, {0,1}};
  vector<Point2D> ipoints { {0.25, 0.25}, {0.75, 0.25}, {0.5, 0.5}, {0.25, 0.75}, {0.75, 0.75}};

  vector<Point2D> all_points(bpoints);
  all_points.insert(all_points.end(), ipoints.begin(), ipoints.end());

  const double vdist = 0.9;
  vector<Triangle> tris;
  // vector<Triangle> tris = triangulate_domain(&all_points[0], (uint)bpoints.size(),
  // 					     (uint)all_points.size(), vdist);

  // cout << "Points: " << endl;
  // for (const auto p : all_points)
  //   cout << p[0] << ", " << p[1] << '\n';
  // cout << endl;
  // cout << "Tris: " << endl;
  // for (const auto t : tris)
  //   cout << t[0] << " " << t[1] << " " << t[2] << '\n';
  // cout << endl << endl;


  // non-convex case
  bpoints = vector<Point2D> { {0,0}, {2,0}, {2,2}, {1,2}, {1,1}, {0,1} };
  ipoints = vector<Point2D> { {0.5, 0.5}, {1, 0.5}, {1.5, 0.5}, {1.5, 1}, {1.5, 1.5}};
  all_points = bpoints;
  all_points.insert(all_points.end(), ipoints.begin(), ipoints.end());

  tris = triangulate_domain(&all_points[0], (uint)bpoints.size(),
			    (uint) all_points.size(), vdist);

  cout << "--- Non-convex case: ---" << endl<<endl;
  cout << "Points: " << endl;
  for (const auto p : all_points)
    cout << p[0] << ", " << p[1] << '\n';
  cout << endl;
  cout << "Tris: " << endl;
  for (const auto t : tris)
    cout << t[0] << " " << t[1] << " " << t[2] << '\n';
  cout << endl << endl;
}

// ----------------------------------------------------------------------------
void test_fit_to_plane()
// ----------------------------------------------------------------------------  
{
  vector<Point3D> pts { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  //vector<Point3D>  pts { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 0.01}};
  
  const auto plane_eq = fit_points_to_plane(&pts[0], (uint)pts.size());

  cout << "Plane equation is: " << endl;
  cout << plane_eq[0] << " x + " << plane_eq[1] << " y + " << plane_eq[2] << " z + ";
  cout << plane_eq[3] << " = 0" << endl;
  
  vector<double> dists(pts.size());
  transform(pts.begin(), pts.end(), dists.begin(), [&](const Point3D& p) {
      return p[0] * plane_eq[0] + p[1] * plane_eq[1] + p[2] * plane_eq[2] + plane_eq[3];
    });
  cout << "Distances to plane are: " << endl;
  copy(dists.begin(), dists.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_volume_tesselation()
// ----------------------------------------------------------------------------
{
  // Tesselate regular prism
  // const vector<Point3D> corners { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, 
  //                               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1} };
  // const vector<Segment> edges   { {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom face
  //                               {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top face
  //                               {0, 4}, {1, 5}, {2, 6}, {3, 7}};  // "vertical" edges
  // const vector<FaceLoop> faces { { {0, 1, 2, 3}, false}, // bottom (zmin) face
  //                                { {4, 5, 6, 7}, true}, // top (zmax) face
  //                                { {0, 9, 4, 8}, true}, // front (ymin) face
  //                                { {2, 11, 6, 10}, true}, // back (ymax) face
  //                                { {8, 7, 11, 3}, true}, // left (xmin) face
  //                                { {1, 10, 5, 9}, true}}; // right (xmax) face

  
  // Tesselate concave object
  const vector<Point3D> corners { {0, 0, 0}, {1, 0, 0}, {1, 1, 0},
                                  {0.5, 1, 0}, {0.5, 0.5, 0}, {0, 0.5, 0},
                                  {0, 0, 1}, {1, 0, 1}, {1, 1, 1},
                                  {0.5, 1, 1}, {0.5, 0.5, 1}, {0, 0.5, 1} };

  const vector<Segment> edges { {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
                                {6, 7}, {7, 8}, {8, 9}, {9, 10}, {10, 11}, {11, 6},
                                {0, 6}, {1, 7}, {2, 8}, {3, 9}, {4, 10}, {5, 11} };

  const vector<FaceLoop> faces { { {0, 1, 2, 3, 4, 5}, false}, // bottom
                                 { {6, 7, 8, 9, 10, 11}, true}, // top
                                   { { 0, 13, 6, 12}, true}, // front
                                   { { 1, 14, 7, 13}, true}, // right
                                   { { 2, 15, 8, 14}, true}, // back 1
                                   { { 3, 16, 9, 15}, true}, // left 1
                                   { { 4, 17, 10, 16}, true}, // back 2
                                   { { 5, 12, 11, 17}, true} };
  
  // const vector<Point3D> corners { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0} };
  // const vector<Segment> edges   { {0, 1}, {1, 2}, {2, 3}, {3, 0}};
  // const vector<FaceLoop> faces { { {0, 1, 2, 3}, false}};
  
  const SimpleVolumeType v {};
  SimplePolyhedron spoly {corners, edges, faces, v};

  cout << "Polyhedron:\n" << endl;
  cout << spoly << endl << endl;

  const double vdist = 0.2; //0.3; //3; //0.5; //0.07; //0.2;
  spoly.tesselate(vdist);

  cout << "Writing wirefame: " << endl;
  spoly.writeTesselatedOutline(cout);

  cout << "Writing shell: " << endl;
  spoly.writeTesselatedShell(cout);

  cout << "Writing tets: " << endl;
  spoly.writeTesselatedVolume(cout);

}

// ----------------------------------------------------------------------------
void test_linear_system()
// ----------------------------------------------------------------------------
{
  const vector<double> m {1, -2, 0, 0, 8, 6, -7,  -3, 1};
  const vector<double> rhs {1, 2, 3};
  vector<double> result(rhs.size());

  // using LU
  solve_linear_system<3>(&m[0], &rhs[0], &result[0]);

  copy(result.begin(), result.end(), ostream_iterator<double>(cout, "\n"));
  cout << endl << endl;
  // using Cramer's rule
  const Point3D c1 {m[0], m[1], m[2]};
  const Point3D c2 {m[3], m[4], m[5]};
  const Point3D c3 {m[6], m[7], m[8]};
  const Point3D rhs2 {1, 2, 3};
  
  const Point3D res = solve_3D_matrix(c1, c2, c3, rhs2);
  cout << res;

  // timing issues
  const int N = 1000000;

  cout << "Now timing " << N << " separate solutions of the system." << endl;
  CPUclock clock;
  clock.initTime();

  for (uint i = 0; i != N; ++i)
    solve_3D_matrix(c1, c2, c3, rhs2);

  cout << clock.getInterval() << " seconds elapsed with Cramer. " << endl;
  
  clock.initTime();
  for (uint i = 0; i != N; ++i)
    solve_linear_system<3>(&m[0], &rhs[0], &result[0]);

  cout << clock.getInterval() << " seconds elapsed with LU." << endl;
}

// ----------------------------------------------------------------------------
void test_polyhedron3D_tesselation()
// ----------------------------------------------------------------------------
{
  // cube corners
  // const vector<Point3D> bpoints { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
  //                                 {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  
  // const vector<Triangle> btris {
  //   {0, 3, 1}, {0, 2, 3}, // bottom
  //   {4, 5, 7}, {4, 7, 6}, // top
  //   {0, 4, 2}, {4, 6, 2}, // left
  //   {1, 3, 5}, {3, 7, 5}, // right
  //   {0, 5, 4}, {0, 1, 5}, // front
  //   {2, 7, 3}, {2, 6, 7}}; // back
  // const double vdist = 0.25; //0.5;//1;//0.5;

  // simplex corners
  // const vector<Point3D> bpoints { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  // const vector<Triangle> btris { {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3} };
  // const  double vdist = 0.2;//0.15;//0.15;
  
  // // extruded "hexagon"
  const double C = sqrt(3)* 0.5;
  const vector<Point3D> bpoints {
    {1, 0, 0}, {0.5, C, 0}, {-0.5, C, 0}, {-1, 0, 0}, {-0.5, -C, 0}, {0.5, -C, 0},
    {1, 0, 1}, {0.5, C, 1}, {-0.5, C, 1}, {-1, 0, 1}, {-0.5, -C, 1}, {0.5, -C, 1}};
  const vector<Triangle> btris 
    { {0, 6, 5}, {6, 11, 5},
      {1, 7, 0}, {7,  6, 0},
      {2, 8, 1}, {8, 7, 1},
      {2, 9, 8}, {9, 2, 3},
      {3, 4, 10}, {3, 10, 9},
      {4, 5, 10}, {5, 11, 10},
      {11, 6, 7}, {11,7,8}, {11, 8, 10}, {10, 8, 9},
      {5, 1, 0}, {5, 2, 1}, {5, 4, 2}, {4, 3, 2}};
    const double vdist = 0.2;
                                
  const auto mesh = tesselatePolyhedron3D(&bpoints[0], (uint)bpoints.size(),
                                          &btris[0], (uint)btris.size(),
                                          vdist);
  cout << "Number of mesh elements: " << mesh.tets.size() << endl;

  ofstream os_points("mesh_points.mat");
  ofstream os_tets("mesh_tets.mat");
  for (const auto p : mesh.points) os_points << p;
  for (const auto t : mesh.tets) os_tets << t;
  os_points.close();
  os_tets.close();
  
}
// ----------------------------------------------------------------------------
void test_circumsphere()
// ----------------------------------------------------------------------------
{
  // base triangle
  array<Point3D, 3> tri { Point3D { 0, 0, 0}, Point3D { 1, 0, 0}, Point3D{1, 1, 0} };

  // bunch of random points
  const int N = 3000;
  vector<double> coords(3 * N);
  generate(coords.begin(), coords.end(), [] () {return random_uniform(-500, 500);});
  vector<Point3D> points;
  for (uint i = 0; i != N; ++i) 
    points.push_back( {coords[3*i], coords[3*i+1], coords[3*i+2]});

  // testing circumsphere on the generated points
  Point3D center;
  double rad2;
  for (uint i = 0; i != N; ++i) 
    if (fitting_sphere(tri[0], tri[1], tri[2], points[i], center, rad2)) {
      if (!(center[2] * points[i][2] > 0)) {
        cout << "Assumption failed for point " << i << ", which is: " << points[i];
        cout << " and has center: " << center << endl;
      }
    } else {
      throw runtime_error("failed to fit sphere");
    }
}

// ----------------------------------------------------------------------------
void test_triangle_intersection()
// ----------------------------------------------------------------------------
{
  // enum IsectCase {
  //   DISJOINT,         // no intersection
  //   AT_POINT,         // share a common point only
  //   ALONG_BOUNDARY,   // boundary of one intersects interior or boundary of other
  //   OVERLAPPING,      // interiors of A and B intersect
  //   CONTAINING,       // one fully contained in other
  //   COLINEAR_OVERLAP  // two 2D segments overlap colinearly 
  // };

  
  vector<Point3D> ref = { {0,0,0}, {1, 0, 0}, {1, 1, 0}};
  vector<Point3D> t1 = { {0,0,0}, {0.9,0,0}, {1, 1.1, 0}}; // should be overlapping (3)
  vector<Point3D> t2 = { {0,0,0.001}, {1, 0, 0.001}, {1, 1, 0.001}};// should not intersect (0)
  vector<Point3D> t3 = { {0,0,0.01}, {1, 0, -0.01}, {1, 1, 0.01}}; // should cut across (3)
  vector<Point3D> t4 = { {0, 0, 0}, {1.1, 1.1, 0}, {0, 1, 0}}; // should be boundary intersect (2)
  vector<Point3D> t5 = { {0.01, 0.01, 0}, {0.99, 0.01, 0}, {0.99, 0.99, 0}}; // contained (4)
  vector<Point3D> t6 = { {0.5, 0.5, 0}, {1, 0, 0}, {1, 0, 1}}; // should cut  across (2)
  vector<Point3D> t7 = { {0.5, 0, -1}, {1, 1, 0}, {0.5, 0, 1}}; // should overlap (2)

  const double tol = 1e-6;
  cout << "t1: " << isect_triangle_triangle_3D(&ref[0], &t1[0], tol) << endl;
  cout << "t2: " << isect_triangle_triangle_3D(&ref[0], &t2[0], tol) << endl;
  cout << "t3: " << isect_triangle_triangle_3D(&ref[0], &t3[0], tol) << endl;
  cout << "t4: " << isect_triangle_triangle_3D(&ref[0], &t4[0], tol) << endl;
  cout << "t5: " << isect_triangle_triangle_3D(&ref[0], &t5[0], tol) << endl;
  cout << "t6: " << isect_triangle_triangle_3D(&ref[0], &t6[0], tol) << endl;
  cout << "t7: " << isect_triangle_triangle_3D(&ref[0], &t7[0], tol) << endl;
  
}

// ----------------------------------------------------------------------------
void test_solve_linear_system()
// ----------------------------------------------------------------------------
{
  const int N = 2000;
  vector<double> m(N*N, 0);
  vector<double> rhs(N,0);
  vector<double> result(N,0);
  generate(m.begin(), m.end(), []() {return random_uniform(-1, 1);});
  generate(rhs.begin(), rhs.end(), []() {return random_uniform(-1, 1);});

  const bool succeed = solve_linear_system<N>(&m[0], &rhs[0], &result[0]);
}

// ----------------------------------------------------------------------------
void test_clip_grid_2D()
// ----------------------------------------------------------------------------
{
  // vector<Point2D> corners = {{0, 1}, {1, 0}, {2, 1}, {1, 2}};
  // vector<Point2D> corners = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  vector<Point2D> corners = {{0,1}, {0.5, 1}, {0.5, 0}, {1, 0}, {1, 1}, {1.5, 1},
                             {0.75, 2}};
  
  const double vdist = 0.21;
  const uint res_x = 43;
  const uint res_y = 43;
  
  const ClippedGrid<2> cg = clip_grid_polygon_2D(&corners[0], (uint)corners.size(),
                                                 vdist, res_x, res_y);

  cout << "bbox:\n";
  copy(cg.bbox.begin(), cg.bbox.end(), ostream_iterator<double>(cout, ", "));
  cout << "\nres:\n";
  copy(cg.res.begin(), cg.res.end(), ostream_iterator<double>(cout, ", "));
  cout << "\ntype:\n";
  copy(cg.type.begin(), cg.type.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_clip_grid_3D()
// ----------------------------------------------------------------------------
{
  //vector<Point3D> corners = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0} };
  vector<Point3D> corners = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  Triangle tri = {0, 1, 2};

  uint res_x = 20;//70; //4;
  uint res_y = 20;//70; //4;
  uint res_z = 20;//70; //8;
  
  const auto cg = clip_grid_shell_3D(&corners[0], 3, &tri, 1, 0.1, res_x, res_y, res_z);

  cout << "bbox:\n";
  copy(cg.bbox.begin(), cg.bbox.end(), ostream_iterator<double>(cout, ", "));
  cout << "\nres:\n";
  copy(cg.res.begin(), cg.res.end(), ostream_iterator<double>(cout, ", "));
  cout << "\ntype:\n";
  copy(cg.type.begin(), cg.type.end(), ostream_iterator<double>(cout, ", "));
  cout << endl;  
}

// ----------------------------------------------------------------------------
void test_parametric_volume_tesselation()
// ----------------------------------------------------------------------------
{
  //string filename = "data/trimmed_cube_tri.g22";
  string filename = "data/pct12108_fem_1couple_model_0_simplified_trim3.g22";
  VolumeModelFileHandler filehandler;

  shared_ptr<ftVolume> testVolume = filehandler.readVolume(filename.c_str());

  GoParametricTesselableVolume ptvolume(*testVolume);

  //const double vdist = 0.007; //0.007; //0.015; //0.02; //0.01; //0.005;
  //const double vdist = 0.02; // gives trouble generating tets
  //const double vdist = 0.04; // OK
  //const double vdist = 0.03; // OK
  //const double vdist = 0.003; // OK  - 100.000 cells
  //const double vdist = 0.01; // for profiling

  //const double vdist = 10; // for second model: (single) ipoint ends up outside shell

  //const double vdist = 1.3;
  //const double vdist = 3.5;

  //const double vdist = 0.9; // make it work for this value
  //const double vdist = 0.7;
  //const double vdist = 1.5;

  //const  double vdist = 2;
  //const double vdist = 1;
  const double vdist = 0.5;
  //const double vdist = 1; // @@ leads to trouble with non-delaunay tets
  //const double vdist = 0.3;// @@ leads to trouble with surface triangulation
  //const double vdist = 0.4;
  
  
  
  ptvolume.tesselate(vdist);

  //ptvolume.writeTesselatedOutline(cout);
  ptvolume.writeTesselatedShell(cout);

  cout << endl;
  cout << "Number of tets generated: " << ptvolume.numTets() << endl;
  
  ptvolume.writeTesselatedVolume(cout);

  vector<GoParametricTesselableVolume::PointType> pts = ptvolume.volumePoints();
  vector<Point3D> points_3D;
  for (auto p : pts)
    points_3D.push_back(Point3D {p.pos[0], p.pos[1], p.pos[2]});
  
  vector<int> non_del_tets = identify_nondelaunay_tets(points_3D,
                                                       ptvolume.getTets(),
                                                       ptvolume.numFacePoints());

  uint count = 0;
  for (size_t i = 0; i != non_del_tets.size(); ++i)
    if (non_del_tets[i] >= 0)
      ++count;
  
  cout << "number of non-delaunay tets: " << count << endl;
  vector<uint> int_tets = interior_tets(ptvolume.getTets(), ptvolume.numFacePoints());
  cout << "number of interior tets: " << accumulate(int_tets.begin(), int_tets.end(), 1) << endl;

  vector<pair<uint, uint>> int_nondel_tets;
  
  for (size_t i = 0; i != non_del_tets.size(); ++i)
    if ( (non_del_tets[i] >= 0) && (int_tets[i] == 1) )
      int_nondel_tets.push_back( {i, non_del_tets[i]} );
  
  cout << "Non-delaunay interior tets: " << endl;
  for (size_t i = 0; i != int_nondel_tets.size(); ++i)
    cout << "Tet: " << int_nondel_tets[i].first << " has point: " << int_nondel_tets[i].second << endl;
}

};
