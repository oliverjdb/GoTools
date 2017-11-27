#ifndef _TESSELATE_UTILS_H
#define _TESSELATE_UTILS_H

#include <array>
#include <vector>
#include <algorithm>
#include "common_defs.h"

namespace TesselateUtils {

const double PI = 3.14159265358979323;

// ----------------------------------------------------------------------------  
template<typename P>
inline double dist2(const P& p1, const P& p2)
// ----------------------------------------------------------------------------
{
  return p1.dist2(p2);
}

// ----------------------------------------------------------------------------
template<typename P>
inline double dist(const P& p1, const P& p2) 
// ----------------------------------------------------------------------------
{
  return std::sqrt(dist2(p1, p2));
}

// ----------------------------------------------------------------------------
template<typename P>
inline void cross(const P& u, const P& v, P& result)
// ----------------------------------------------------------------------------
{
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}
  
// ----------------------------------------------------------------------------
template<typename P>
std::vector<P> interpolate(const P& p1, const P& p2, unsigned int num,
                           bool include_endpoints = true);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns signed area of polygon.  A counterclockwise-ordered polygon has
// positive area.
template<typename P>
double polygon_area(const P* const poly, const unsigned int num_corners);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------  
template<int Dim>
std::array<double, 2*Dim> compute_bounding_box(const PointXD<Dim>* const points,
                                               const unsigned int num_points);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Returns a bounding box on the form [xmin, xmax]
template<typename P>  
std::array<double, 2> bounding_box_1D(const P* const points,
                                      const unsigned int num_points);
// ----------------------------------------------------------------------------  
  
// ----------------------------------------------------------------------------
// Returns a bounding box on the form [xmin, xmax, ymin, ymax]
template<typename P>  
std::array<double, 4> bounding_box_2D(const P* const points,
                                      const unsigned int num_points);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns a bounding box on the form [xmin, xmax, ymin, ymax, zmin, zmax]
template<typename P>  
std::array<double, 6> bounding_box_3D(const P* const points,
                                      const unsigned int num_points);
// ----------------------------------------------------------------------------  
  
// ----------------------------------------------------------------------------
// Generate grid of points with resolution (nx+2) x (ny+2) in the rectangle
// defined by the two opposite corners c1 and c2.  nx and ny denote the number
// of _internal_ grid points, hence the (nx+2) x (ny+2) resolution
template<typename P>
std::vector<P> generate_grid_2D(const P& c1,
                                const P& c2,
                                unsigned int nx,
                                unsigned int ny,
                                bool include_boundary_points = true);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Generate grid of points with resolution (nx+2) x (ny+2) x (nz+2) in the prism
// defined by the two opposite corners c1 and c2.  nx, ny and nz denote the
// number of _internal_ grid points, hence the (nx+2) x (ny+2) x (nz+2)
// resolution
template<typename P>
std::vector<P> generate_grid_3D(const P& c1,
                                const P& c2,
                                unsigned int nx,
                                unsigned int ny,
                                unsigned int nz,
                                bool include_boundary_points = true);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------  
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& arg);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Check if a point is inside a polygon
template<typename P> bool inpolygon(const P& pt, const P* const poly, 
				    const unsigned int num_corners,
                                    const double tol,
                                    bool& on_boundary);
// ----------------------------------------------------------------------------      
  
// ----------------------------------------------------------------------------
// Return the points that are inside the polygon
template<typename P>
std::vector<P> inpolygon(const P* const pts, const unsigned int num_pts,
			 const P* const poly, const unsigned int num_corners, 
			 const double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Compute approximate normal of a polygon in 3D space
template<typename P>
P compute_polygon_normal(const P* const pts, const unsigned int num_pts);
// ----------------------------------------------------------------------------
  
  
// ----------------------------------------------------------------------------
// Check if a 3D point is inside a closed shell consisting of triangles
template<typename P, typename T>
bool inside_shell(const P& pt, const P* const bpoints, const T* const tris,
                  const unsigned int num_tris, const double tol, bool& on_bnd);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Returns the 3D points that are inside the closed shell defined by a set of
// triangles
template<typename P, typename T>
std::vector<P> inside_shell(const P* const pts, const unsigned int num_pts,
                            const P* const bpoints, const T* const tris,
                            const unsigned int num_tris, const double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Find the point on the triangulated surface that is closest to the point 'p'.
// The index to one triangle containing this point (which may not be unique) is
// returned in 'tri_ix'.
template<typename P, typename T>
P closest_point_on_triangle_surface(const P& p,
                                    const P* const surf_pts,
                                    const unsigned int num_surf_pts,
                                    const T* const surf_tris,
                                    const uint num_tris,
                                    uint& tri_ix);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------  
template<typename P>
P closest_point_on_loop(const P& p,
                        const P* const bpoints,
                        const unsigned int num_bpoints,
                        uint& seg_ix);
// ----------------------------------------------------------------------------  

  
// ----------------------------------------------------------------------------
// 'p' is the point we want to compute distance from.  'a' and 'b' are two
// distinct points on the infinite line (thus defining the line).  The returned
// distance is signed: a positive number indicates that the point is on the
// "right" of the directed line from 'a' to 'b'.
template<typename P> double projected_distance_to_line_2D(P p, P a, P b);
template<typename P> double projected_distance_to_line_3D(P p, P a, P b);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Return the projected distance from the point to the plane defined by the
// three points in 'tricorners'.  The normalized direction vector from the point
// to the plane will be written to 'dir'.
template<typename P> inline
double projected_distance_to_plane(const P& pt, const P* const tricorners, P& dir);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------    
// Check if the point p is within distance 'tol' of the line segment defined by
// points a and b.  'in_3D' should be 'true' for 3D, and 'false' for 2D
template<typename P> inline
bool point_on_line_segment(const P& p, const P& a, const P& b, double tol, bool in_3D);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------    
// Check if a point is within distance 'tol' from a triangle in 3D space.
template<typename P> inline
bool point_on_triangle(const P& p, const P& c1, const P& c2, const P& c3,
                       const double tol);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template<typename P> inline
P triangle_normal(const P* const tripoints, double& area);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
enum ProjectionType { // whether the projected point on the triangle falls
  INSIDE,             // inside the triangle, on its edge or in one of its corners
  ON_EDGE,
  AT_CORNER
};
// Report whether or not a given point 'p' is within distance 'tol' of a
// triangle specified by corners 'tripoints'.  If it is the case, a certain
// amount of additional information is provided.  This includes the position of
// the projected point itself ('proj_pt'), the signed, normal distance from 'p' the
// plane containing the triangle, the projection type (whether 'proj_pt' is
// found inside the triangle, along one of its edges or at one of its corners),
// and in the case of a projected point located on an edge, the normalized
// direction vector of that edge.
template<typename P> inline
bool point_on_triangle(const P& p, const P* const tripoints,
                       const double tol, P& proj_pt, double& dist_to_plane,
                       ProjectionType& ptype, P& edge_dir);
// ----------------------------------------------------------------------------      

// ----------------------------------------------------------------------------      
template<typename P> inline
bool triangle_area_3D(const P& c1, const P& c2, const P& c3);
// ----------------------------------------------------------------------------      


// ----------------------------------------------------------------------------
// Check if a perpendicular projection of a 3D point onto the plane of the
// triangle falls inside the triangle.  If so, the output variables 'dist_2' and
// 'sign' will respectively inform about the squared distance from the point to
// the triangle, and whether the point is on the 'back' (positive sign), or
// 'front' (negative sign) of the oriented triangle.
template<typename P> inline
bool projects_to_triangle(const P& pt, const P* const tricorners, const double tol,
                          double& dist_2, int& sign);
// ----------------------------------------------------------------------------      
  
// ----------------------------------------------------------------------------
// Returns 'true' if the line through point 'p' that perpendicularly intersects
// the line passing through 'a' and 'b' has its itnersection point within the
// segment ab.
template<typename P> inline
bool projects_to_segment(const P& p, const P& a, const P&b);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------    
template<typename P, int Dim> inline  
P projected_point_on_segment(const P& p, const P& a, const P&b, bool& at_corner);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// This function is only implemented for (and only makes sense for) 3D points.
// It checks whether the ray refined by the point 'pt' and direction 'dir'
// intersects with the triangle defined by the points 'p1, 'p2' and 'p3.  If so,
// the function returns true.  Moreover, the squared distance to the
// intersection is returned in the 'dist_2' output argument, and the sign of the
// distance is returned in 'sign'.  A distance is positive if the ray and the
// triangle normal point in the same direction (positive scalar product).  The
// orientation of the triangle's normal is defined by the order of the triangle
// corners, which should make a counterclockwise loop.  If the ray is parallel
// to the plane of the triangle, it is interpreted not to intersect, regardless
// of whether it is in the plane or not.  However, this case will be flagged by
// setting 'sign' to zero, whereas if the non-intersection was caused by the
// intersection point falling outside the triangle, the 'sign' will be nonzero
// (but otherwise meaningless).
template<typename P> inline
bool ray_intersects_face(const P& pt, const P& dir,
                         const P& p1, const P& p2, const P& p3,
                         double tol,
                         double& dist_2, int& sign);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Check if 3D point 'pt' is on the "inside" or the "outside of the oriented
// plane defined by the oriented triangle face given by its corner points in
// counterclockwise order, 'p1', 'p2', and 'p3'.
template<typename P>
bool point_on_inside_of_face(const P& pt, const P& p1, const P& p2, const P& p3,
                             const double tol);
// ----------------------------------------------------------------------------  
  
// ----------------------------------------------------------------------------
// return 'true' if the angle abc is acute
template<typename P> inline
bool acute_angle(const P& a, const P& b, const P& c);
// ----------------------------------------------------------------------------      

// ----------------------------------------------------------------------------
// Mirror point 'p' accross the line passing through 'a' and 'b', and return the
// result.
template<typename P> inline P mirror_point_2D(P p, const P& a, const P& b);
// ----------------------------------------------------------------------------

  
// ----------------------------------------------------------------------------
// mirror the provided set of points across the line passing through a and b,
// and return the result (the mirror points) as a vector
template<typename P> inline
std::vector<P> mirror_points_2D(const P* const pts,
				const unsigned int num_pts,
				const P& a,
				const P& b);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Mirror point 'p' across the plane containing the triangle whose corners are
// given in 'tricorners', and return the result.
template<typename P> inline P mirror_point_3D(P p, const P* const tricorners);
// ----------------------------------------------------------------------------
  
// ----------------------------------------------------------------------------
// mirror the provided set of points across the triangle given by the three
// corners pointed to in 'tricorners', and pass the result (the mirror points)
// as a vector.
template<typename P> inline
std::vector<P> mirror_points_3D(const P* const pts,
                                const unsigned int num_pts,
                                const P* const tricorners);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------    
double random_uniform(double minval, double maxval);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// similar to MATLAB's "find".  There are two versions, depending on whether
// clarity of code (result in return value) or speed (result as output argument)
// is favored.
template<typename T> inline
std::vector<unsigned int> locate_nonzeros(const std::vector<T> vec);

template<typename T> inline
void locate_nonzeros(const std::vector<T> vec, std::vector<unsigned int>& target);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Extract all elements from a range that satisfy a predicate. The predicate
// should take a reference to T and return a bool (or int).  The first vector in
// the result contains the extracted values, the second vector contain the
// respective indices.
template<typename T, typename Pred> inline
std::pair<std::vector<T>, std::vector<unsigned int>>
extract_from_range(const T* const range_start,
		   unsigned int range_length, const Pred& fun);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
template<typename T>
std::vector<T> extract_from_range(const std::vector<T>& vec,
				  const std::vector<unsigned int>& ixs)
// ----------------------------------------------------------------------------
{
  std::vector<T> res(ixs.size());
  transform(ixs.begin(), ixs.end(), res.begin(), [&vec](unsigned int i) {return vec[i];});
  return res;
}
  
// ----------------------------------------------------------------------------
// solve a 2x2 linear system.  The columns of the 2x2 matrix are given in the
// first two arguments, the right-hand-side in the third.  If matrix is
// singular, the result will contain NaN-values.
template<typename P> inline
P solve_2D_matrix(const P& Mcol1, const P& Mcol2, const P& rhs);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// solve a 3x3 linear system.  The columns of the 3x3 matrix are given in the
// first three arguments, the right-hand-side in the fourth.  If matrix is
// singular, the result will contain NaN-values.
template<typename P> inline
P solve_3D_matrix(const P& Mcol1, const P& Mcol2, const P& Mcol3, const P& rhs);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Solve N-sized linear system.  Result vector will be written to array pointed
// to by 'result'.  Function returns 'true' if success.
template<int N> inline
bool solve_linear_system(const double* const m,
                         const double* const rhs,
                         double* const result);
// ----------------------------------------------------------------------------  

  
// ----------------------------------------------------------------------------
// Compute 3D determinant of matrix defined by the columns Mcol1, Mcol2 and
// Mcol3.
template<typename P> inline
double determinant3D(const P& Mcol1, const P& Mcol2, const P& Mcol3);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Compute the center point and the radius of the circle that circumscribes the
// triangle defined by p1, p2 and p3.  If the points are collinear, return
// 'false'.  The circle center and squared radius are returned in the two last
// (output) arguments.
template<typename P> inline
bool circumscribe_triangle(const P& p1, const P& p2, const P& p3,
			   P& circ_center, double& radius2);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// Determine if the triangle defined by p1, p2 and p3 is obtuse
// (largest angle > 90 degrees).
template<typename P> inline
bool is_obtuse_triangle(const P& p1, const P& p2, const P& p3);
// ----------------------------------------------------------------------------
  
// ----------------------------------------------------------------------------
// Compute the unique sphere that contains the four provided points on its
// boundary.  Return false if no such sphere exist (i.e. points are planar)
template<typename P> inline
bool fitting_sphere(const P& p1, const P& p2, const P& p3, const P& p4,
                    P& center, double& radius2);
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// check whether two segments intersect.  A positive value for the tolerance
// means that segments within a distance of 'tol' times their length scale
// counts as an intersection. A negative value for the tolerance, on the other
// hand, means that each line segment must cross the other by a length of at
// least one times |tol| times their length scale.  In this case, a mere touch
// is not sufficient to count as a intersection.
template<typename P> inline
bool segments_intersect_2D(const P& seg1_a, const P& seg1_b,
			   const P& seg2_a, const P& seg2_b, const double tol);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Check whether two triangles (p1, p2, p3) and (q1, q2, q3) intersect.  A
// positive value for the tolerance means that triangles closer to each other
// than 'tol' times their length scale are interpreted to intersect.  A negative
// value for the tolerance means there need to be an actual intersection of at
// least their length scale times |tol|.  In this case, a mere touch is not
// sufficient to count as an intersection.
template<typename P> inline
bool triangles_intersect_3D(const P& p1, const P& p2, const P& p3,
                            const P& q1, const P& q2, const P& q3, const double tol);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Check if a line segments intersects a triangle in 3D space.  A positive value
// for the tolerance means that a segment within a distance of 'tol' to the face
// is interpreted to intersect it.  A negative value of 'tol' needs that the
// segment must penetrate the triangle by at least |tol| times the length scale.
template<typename P> inline
bool segment_intersects_face(const P& seg_a, const P& seg_b,
                             const P& tri_a, const P& tri_b, const P& tri_c,
                             const double tol);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Transform a set of 3D points to a set of 2D points, by projecting on the
// plane defined by the first three points in the set, and using the first point
// as the origin in the plane.
template<typename P2D, typename P3D>
std::vector<P2D> transform_to_2D(const std::vector<P3D>& pts);
// ----------------------------------------------------------------------------

}; // end namespace TesselateUtils

#include "tesselate_utils_impl.h"

#endif
