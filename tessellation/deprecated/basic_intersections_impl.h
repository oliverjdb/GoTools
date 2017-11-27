#ifndef _BASIC_INTERSECTIONS_IMPL_H
#define _BASIC_INTERSECTIONS_IMPL_H

#include <vector>
#include <array>
#include "common_defs.h"
#include "tesselate_utils.h"

namespace TesselateUtils {

namespace { // helpers
// ----------------------------------------------------------------------------  
inline bool seg_inside_1D(const double s1_min, const double s1_max,
                          const double s2_min, const double s2_max, const double tol)
// ----------------------------------------------------------------------------  
{
  return ((s1_min > s2_min - tol) && (s1_max < s2_max + tol));
}
  
  
// ----------------------------------------------------------------------------
template<typename P>
inline std::vector<double> reduce_to_1D(const std::vector<P>& points2D)
// ----------------------------------------------------------------------------
{
  // Points are assumed to be collinear. Take first point as origin, and
  // distance to this point as the position along the axis (sign determined by
  // comparing with second point)
  assert(points2D.size() > 1);
  std::vector<double> result(points2D.size(), 0);
  const P origin = points2D[0];
  const P dirvec {points2D[1][0] - origin[0], points2D[1][1] - origin[1]};
  std::transform(points2D.begin(), points2D.end(), result.begin(),
                 [&origin, &dirvec] (const P& p) {
                   const double val = sqrt(dist2(p, origin));
                   const P d = {p[0] - origin[0], p[1] - origin[1]}; 
                   const double sprod = d[0] * dirvec[0] + d[1] * dirvec[1];
                   return (sprod > 0) ? val : -val;
                 });
  
  return result;
}

// ----------------------------------------------------------------------------
template<typename P>
inline std::vector<std::array<double, 2>>
reduce_to_2D(const std::vector<P>& points3D, P& origin, P& v1, P& v2)
// ----------------------------------------------------------------------------
{
  // Points are assumed to be complanar. Take first point as origin, normalized
  // vector from first to second point to be the first axis, and normalized
  // vector, and the second axis to be obtained from the (negative) cross
  // product of the first axis and the plane normal.
  assert(points3D.size() > 2);

  origin = points3D[0];
  v1 = (points3D[1] - points3D[0]) / dist(points3D[1], points3D[0]);
  P tmp = (points3D[2] - points3D[0]);
  P n; cross(v1, tmp, n);
  n /= norm(n);
  v2; cross(n, v1, v2);

  // we now have the orthonormal system (v1, v2, n)
  std::vector<std::array<double,2>> result(points3D.size());
  std::transform(points3D.begin(), points3D.end(), result.begin(),
                 [v1, v2, n, origin] (P p) {
                   p -= origin;
                   P tmp = solve_3D_matrix(v1, v2, n, p);
                   return std::array<double, 2> {tmp[0], tmp[1]};
                 });
  return result;
}
  
// ----------------------------------------------------------------------------
template<typename P>
bool point_inside_triangle(const P& p1, const P* const tri, const double tol)
// ----------------------------------------------------------------------------      
{
  return (projected_distance_to_line_2D(p1, tri[0], tri[1]) < tol) &&
         (projected_distance_to_line_2D(p1, tri[1], tri[2]) < tol) &&
         (projected_distance_to_line_2D(p1, tri[2], tri[0]) < tol);
}

// ----------------------------------------------------------------------------
inline bool bounding_boxes_overlap(const std::array<double, 4>& box1,
                                   const std::array<double, 4>& box2,
                                   const double tol)
// ----------------------------------------------------------------------------
{
  double dummy[2];
  return (isect_seg_seg_1D(&box1[0], &box2[0], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[2], &box2[2], tol, dummy) != DISJOINT);
}

// ----------------------------------------------------------------------------
inline bool bounding_boxes_overlap(const std::array<double, 6>& box1,
                                   const std::array<double, 6>& box2,
                                   const double tol)
// ----------------------------------------------------------------------------
{
  double dummy[2];
  return (isect_seg_seg_1D(&box1[0], &box2[0], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[2], &box2[2], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[4], &box2[4], tol, dummy) != DISJOINT);
}

// ----------------------------------------------------------------------------  
template<typename P>
bool tri_inside_tri(const P* tri1, const P* tri2, const double tol)
// ----------------------------------------------------------------------------
{
  return point_inside_triangle(tri1[0], tri2, tol) &&
         point_inside_triangle(tri1[1], tri2, tol) &&
         point_inside_triangle(tri1[2], tri2, tol);
}

// ----------------------------------------------------------------------------  
template<typename P> std::array<P, 2>
scale_segment_beyond_bounding_box(const P* seg, const std::array<double, 4>& bbox)
// ----------------------------------------------------------------------------
{
  const P dp = seg[1] - seg[0];
  const uint comp = (std::fabs(dp[0]) > std::fabs(dp[1])) ? 0 : 1; // which is the largest component

  const double min_target = bbox[comp];
  const double max_target = bbox[comp + 2];

  std::array<P, 2> result {seg[0], seg[1]};

  int expo = 1;
  
  while ( (std::min(result[0][comp], result[1][comp]) >= min_target) ||
          (std::max(result[0][comp], result[1][comp]) <= max_target)) {

    result[0] = seg[0] - (dp * std::pow(2, expo));
    result[1] = seg[1] + (dp * std::pow(2, expo));

    expo++;
  }
  return result;
}

// ----------------------------------------------------------------------------
enum ParallelState { IDENTICAL, PARALLEL, REGULAR };
template<typename P> ParallelState check_coplanar_state(const P* const tri1,
                                                         const P* const tri2,
                                                         double tol)
// ----------------------------------------------------------------------------
{
  P n1; cross(tri1[1] - tri1[0], tri1[2] - tri1[0], n1);
  n1 /= norm(n1);
  P n2; cross(tri2[1] - tri2[0], tri2[2] - tri2[0], n2);
  n2 /= norm(n2);

  P tmp; cross(n1, n2, tmp);
  if (norm(tmp) < tol) {
    // plane are parallel - but are they identical?  take the distance from one
    // point in each plane, and see if it has a non-neglible component along the
    // normal vector.
    P dp = tri1[0] - tri2[0];
    const double sprod = sqrt(dp[0] * n1[0] +
                              dp[1] * n1[1] +
                              dp[2] * n1[2]);

    return (std::fabs(sprod) < tol) ? IDENTICAL : PARALLEL;
  } 
  return REGULAR;
}


// ----------------------------------------------------------------------------
template<typename P> ParallelState check_colinear_state(const P* const seg1,
                                                        const P* const seg2,
                                                        double tol)
// ----------------------------------------------------------------------------
{
  P u = seg1[1] - seg1[0]; u /= norm(u);
  P v = seg2[1] - seg2[0]; v /= norm(v);
  const double tmp = u[0] * v[1] - u[1] * v[0];

  if (std::fabs(tmp) < tol) {
    // lines defined by segments are parallel - but are they identical?
    P w = seg2[0] - seg1[0];
    const double tmp2 = u[0] * w[1] - u[1] * w[0];

    return (std::fabs(tmp2) < tol) ? IDENTICAL : PARALLEL;
  }
  return REGULAR;
}

  

// ----------------------------------------------------------------------------
template<typename P>
std::array<P, 2> compute_plane_intersection(const P* const tri1,
                                            const P* const tri2)
// ----------------------------------------------------------------------------
{
  P n1; cross(tri1[1] - tri1[0], tri1[2] - tri1[0], n1);
  n1 /= norm(n1);
  P n2; cross(tri2[1] - tri2[0], tri2[2] - tri2[0], n2);
  n2 /= norm(n2);

  P dir; cross(n1, n2, dir); // direction vector

  // We now need to find a specific point on the intersection

  // determine index of largest component
  const auto m_el = std::max_element(&dir[0], &dir[0]+3,
                                [](double a, double b) {
                                  return std::fabs(a) < std::fabs(b);});
  const uint m_ix = uint(m_el - &dir[0]);
  const uint i1 = (m_ix + 1)%3;
  const uint i2 = (m_ix + 2)%3;

  const double d1 = tri1[0][0] * n1[0] + tri1[0][1] * n1[1] + tri1[0][2] * n1[2];
  const double d2 = tri2[0][0] * n2[0] + tri2[0][1] * n2[1] + tri2[0][2] * n2[2];

  const std::array<double, 2> rhs {-d1, -d2};

  // determining matrix column 1 for linear system to solve
  const std::array<double, 2> mcol1 {n1[i1], n2[i1]};
  const std::array<double, 2> mcol2 {n1[i2], n2[i2]};
  std::array<double, 2> res = solve_2D_matrix(mcol1, mcol2, rhs);
  
  P p1;  p1[i1] = res[0]; p1[i2] = res[1]; p1[m_ix] = 0;
  return {p1, p1 + dir};
}

}; // end helper namespace

  

// ----------------------------------------------------------------------------
inline
IsectCase isect_point_seg_1D(const double p, const double* const s, const double tol)
// ----------------------------------------------------------------------------  
{
  const double smin = std::min(s[0], s[1]);
  const double smax = std::max(s[0], s[1]);

  if ((p < smin-tol) || (p > smax + tol)) return DISJOINT; // point outside bounding box
  if ((p < smin+tol) || (p > smax - tol)) return AT_POINT; // point "barely inside"

  return CONTAINING; // point well inside
  
}


// ----------------------------------------------------------------------------
inline  
IsectCase isect_seg_seg_1D(const double* const s1,
                            const double* const s2,
                            const double tol,
                            double* const result)
// ----------------------------------------------------------------------------  
{
  const double s1_min = std::min(s1[0], s1[1]);
  const double s1_max = std::max(s1[0], s1[1]);
  const double s2_min = std::min(s2[0], s2[1]);
  const double s2_max = std::max(s2[0], s2[1]);

  if ((s1_max < s2_min - tol) || (s2_max < s1_min - tol)) return DISJOINT; 

  if (s1_max < s2_min + tol) { // s1_max touches s2_min
    result[0] = s2_min; result[1] = s1_max;
    return AT_POINT;
  } else if (s2_max < s1_min + tol) { // s2_max touches s1_min
    result[0] = s1_min; result[1] = s2_max;
    return AT_POINT;
  } else if (seg_inside_1D(s1_min, s1_max, s2_min, s2_max, tol)) {// s1 in side s2
    result[0] = s1_min; result[1] = s1_max; 
    return CONTAINING;
  } else if (seg_inside_1D(s2_min, s2_max, s1_min, s1_max, tol)) {
    result[0] = s2_min; result[1] = s2_max;
    return CONTAINING;
  }
  // if we got here, segments are nontrivially overlapping
  result[0] = std::max(s1_min, s2_min);
  result[1] = std::min(s1_max, s2_max);
  return OVERLAPPING;
}
  
// ----------------------------------------------------------------------------
template <typename P>
IsectCase isect_seg_seg_2D(const P* const s1,
                           const P* const s2,
                           const double tol)
                            
// ----------------------------------------------------------------------------    
{
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_2D(s1, 1), bounding_box_2D(s2, 1), tol))
    return DISJOINT;

  // check if lines are colinear/parallel
  const auto clstate = check_colinear_state(s1, s2, tol);
  if (clstate == IDENTICAL) { 
    // segments lie on the same line.  Reduce dimensionality of problem
    const auto segs_1D = reduce_to_1D(std::vector<P> {s1[0], s1[1], s2[0], s2[1]});
    double res_1D[2];
    const auto icase = isect_seg_seg_1D(&segs_1D[0], &segs_1D[2], tol, res_1D);
    
    return (icase == OVERLAPPING) ? COLINEAR_OVERLAP : icase;
  } else if (clstate == PARALLEL) {
    // lines are parallel but not identical.  Cannot intersect
    return DISJOINT;
  }

  // Segments are not collinear. Compute unique intersection point of the two
  // lines on which segments lie.
  // We determine u and v such that: (1-u) * s1[0] + u * s1[1] = (1-v) *
  // s2[0] + v * s2[1]
  const auto uv = solve_2D_matrix( P {s2[1][0] - s2[0][0], s2[1][1] - s2[0][1]},
                                   P {s1[0][0] - s1[1][0], s1[0][1] - s1[1][1]},
                                   P {s1[0][0] - s2[0][0], s1[0][1] - s2[0][1]});

  const P result = {(s1[0][0] * (1 - uv[0]) + s1[1][0] * uv[0]),  // this is the
                    (s1[0][1] * (1 - uv[0]) + s1[1][1] * uv[0])}; // intersection pt.

  const auto pts_1D = reduce_to_1D( std::vector<P> {result, s1[0], s1[1], s2[0], s2[1]} );

  const auto icase1 = isect_point_seg_1D(pts_1D[0], &pts_1D[1], tol);
  const auto icase2 = isect_point_seg_1D(pts_1D[0], &pts_1D[3], tol);
  
  if ( (icase1 == DISJOINT) || (icase2 == DISJOINT))
    return DISJOINT; // intersection point clearly outside at least one segment

  if ( (icase1 == AT_POINT) || (icase2 == AT_POINT))
    return AT_POINT; // at least one of the segments intersect in an enpoint

  return OVERLAPPING; // if we got here, there was one intersection point that
                      // was clearly inside both segments
}
  
// ----------------------------------------------------------------------------    
template<typename P>
IsectCase isect_seg_triangle_2D(const P* const s, // pointer to two points
                                 const P* const tri, // pointer to three points
                                 const double tol)
// ----------------------------------------------------------------------------  
{
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_2D(s, 2), bounding_box_2D(tri, 3), tol)) 
    return DISJOINT;

  // Check if segment is fully contained within triangle
  if (point_inside_triangle(s[0], tri, tol) && point_inside_triangle(s[1], tri, tol)) {
    return CONTAINING;
  }

  // Segment is not fully contained.  If an intersection exist, there has to be
  // an edge intersection somewhere.
  IsectCase icase[3];
  for (uint i = 0; i != 3; ++i) {
    std::array<P, 2> tseg {tri[i], tri[(i+1)%3]};
    icase[i] = isect_seg_seg_2D(s, &tseg[0], tol);

    if (icase[i] == CONTAINING || icase[i] == COLINEAR_OVERLAP) 
      return ALONG_BOUNDARY;
  }
  // if we got here, there as no colinear overlap between segments.  There may
  // still be AT_POINT or regular OVERLAPs between segments
  
  for (uint i = 0; i != 3; ++i) 
    if (icase[i] == OVERLAPPING)
      return OVERLAPPING;
  
  for (uint i = 0; i != 3; ++i)
    if (icase[i] == AT_POINT)
      return AT_POINT;

  return DISJOINT;
    
}

// ----------------------------------------------------------------------------
template<typename P>
IsectCase isect_line_triangle_2D(const P* const s, // segment defining the line
                                  const P* const tri,
                                  const double tol)
// ----------------------------------------------------------------------------  
{
  // determining a segment large enough that it plays the role of an infinite
  // line in the segment/triangle intersection routine
  const auto longseg = scale_segment_beyond_bounding_box(s, bounding_box_2D(tri, 3));
  return isect_seg_triangle_2D(&longseg[0], tri, tol);
}

// ----------------------------------------------------------------------------  
template<typename P>
IsectCase isect_triangle_triangle_2D(const P* const tri1,
                                      const P* const tri2,
                                      const double tol)
// ----------------------------------------------------------------------------    
{
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_2D(tri1, 3), bounding_box_2D(tri2, 3), tol)) 
    return DISJOINT;
  
  // check if one triangle is fully contained in the other
  if (tri_inside_tri(tri1, tri2, tol) || tri_inside_tri(tri2, tri1, tol))
    return CONTAINING;

  // if we got here, overlap of triangles require overlap of segments.
  // Intersect segments of one triangle against the other triangle
  IsectCase icase[3];
  for (uint i = 0; i != 3; ++i) {
    std::array<P, 2> seg {tri1[i], tri1[(i+1)%3]};
    icase[i] = isect_seg_triangle_2D(&seg[0], tri2, tol);
    assert(icase[i] != CONTAINING);
    if (icase[i] == OVERLAPPING)
      return OVERLAPPING;
  }

  for (uint i = 0; i != 3; ++i) 
    if (icase[i] == ALONG_BOUNDARY)
      return ALONG_BOUNDARY;

  for (uint i = 0; i != 3; ++i) 
    if (icase[i] == AT_POINT)
      return AT_POINT;

  return DISJOINT;
  
}

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, ALONG_BOUNDARY, CONTAINING or OVERLAPPING
template<typename P>
IsectCase isect_triangle_triangle_3D(const P* const tri1,
                                      const P* const tri2,
                                      const double tol)
// ----------------------------------------------------------------------------    
{
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_3D(tri1, 3), bounding_box_3D(tri2, 3), tol))
    return DISJOINT;

  // check if triangles are coplanar
  const auto coplanar_state = check_coplanar_state(tri1, tri2, tol);
  if (coplanar_state == IDENTICAL) {
    P origin, v1, v2; // we will not use these, but function requires them
    const auto system_2D = reduce_to_2D({tri1[0], tri1[1], tri1[2],
                           tri2[0], tri2[1], tri2[2]}, origin, v1, v2);
    
    return isect_triangle_triangle_2D(&system_2D[0], &system_2D[3], tol);
    
  } else if (coplanar_state == PARALLEL) {
    return DISJOINT; // parallel but separate planes cannot intesect
    
  }

  // if we got here, we are in the 'regular' case where triangle planes do intersect
  
  // compute intersection line
  std::array<P, 2> iline = compute_plane_intersection(tri1, tri2);

  const auto icase1 = isect_line_triangle_2D(&iline[0], tri1, tol);
  const auto icase2 = isect_line_triangle_2D(&iline[0], tri2, tol);

  if (icase1 == DISJOINT || icase2 == DISJOINT)
    return DISJOINT;

  if (icase1 == AT_POINT || icase2 == AT_POINT)
    return AT_POINT;

  if (icase1 == ALONG_BOUNDARY || icase2 == ALONG_BOUNDARY)
    return ALONG_BOUNDARY;

  assert(icase1 == OVERLAPPING && icase2 == OVERLAPPING);
  return OVERLAPPING;

}

  
}; // end namespace TesselateUtils

#endif

