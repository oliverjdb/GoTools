#ifndef _TRITOOLS_H
#define _TRITOOLS_H

#include "GoTools/utils/Point.h"
#include <array>
#include <vector>

namespace TriTools {

  // Parametric points specified by 2D Go::Points 
  // Triangles are specified as 3-tuple of ndoe indices.
  struct SurfaceTriangulation
  {
    std::vector<Go::Point> uv;
    std::vector<std::array<int, 3>> triangles;
    size_t num_interior_nodes;
  };
  
  // Return the coefficients [a, b, c] of the implicit equation ax + by + c = 0
  // that describes the line on which the points pt_first and pt_second lie.
  // The coefficients are normalized so that a^2 + b^2 = 1.  With this
  // normalization, for a given point (px, py), the function a * px + b * py + c
  // gives the signed distance from the point to the line.
  std::array<double, 3> normalized_line_equation(const Go::Point& pt_first,
						 const Go::Point& pt_second);

  // Checks which of a number or points is inside a given triangle, for
  // specified margins which may differ from one edge to the next.  The triangle
  // is defined by the three points 'p1', 'p2' and 'p3'.  The points to check
  // against are given in the vector 'points'.  The margins for each edge are
  // specified in 'margin', where 'margin[0]' represents the margin (in terms of
  // Euclidian distance) from the line described by p1->p2, 'margin[1]' the margin
  // for the line described by p2 -> p3 and so on.  The result vector has the
  // same number of entries as 'points', and consists of zeroes and ones.
  // Entries with a one represent points located within the triangle.  The
  // margins should all be positive, and represent the buffer between the line
  // and the interior.
  std::vector<int> points_inside_triangle(const Go::Point& p1,
					  const Go::Point& p2,
					  const Go::Point& p3,
					  const std::vector<Go::Point>& points,
					  const std::array<double, 3>& margin);

  // Identify points located inside the domain bounded by a number of 'loops'
  // (polygons), each specified by a number of points (corners).
  // Counterclockwise loops represent outer boundaries and clockwise loops
  // holes.  The points to test against are provided in the vector 'points', and
  // 'margin' specifies a minimum distance to the boundary (points closer to the
  // boundary than this distance will be considered 'outside').  The result
  // vector has the same number of entries as 'points', and consists of zeroes
  // and ones.  Entries with a one represent points located within the triangle.
  // The margins should all be positive, and represent the buffer between the
  // line and the interior.
  std::vector<int> points_inside_loops(const std::vector<std::vector<Go::Point>>& loops,
				       const std::vector<Go::Point>& points,
				       const double margin);

  // Create a triangulation of a domain bounded by a number of 'loops',
  // (polygons), each specified by a number of points (corners).
  // Counterclockwise loops represent outer boundaries and clockwise loops
  // holes.  The resulting triangulation will be that of a constrained delauney.
  // Each triangle is represented as a triplet of points (first entry in the
  // 'pair' of the result vector components), and a triplet of booleans (second
  // entry in the pair).  The triplet of points represent the triangle corners,
  // and the booleans specify which of the triangle edges are assocated with an
  // external boundary.
  std::vector<std::pair<std::array<Go::Point, 3>, std::array<bool, 3>>>
  triangulate_boundary(const std::vector<std::vector<Go::Point>>& loops);

  SurfaceTriangulation
  triangulate_with_boundaries(const std::vector<Go::Point>& points,
			      const std::vector<std::vector<Go::Point>>& loops);

  
}; // end namespace TriTools

#endif

