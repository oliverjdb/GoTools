#include <vector>
#include <array>
#include "basic_intersections.h"
#include "common_defs.h"
#include "tesselate_utils.h"

using namespace std;
using namespace TesselateUtils;

namespace { // helpers

// ----------------------------------------------------------------------------  
inline bool seg_inside_1D(const double s1_min, const double s1_max,
                          const double s2_min, const double s2_max, const double tol)
// ----------------------------------------------------------------------------  
{
  return ((s1_min > s2_min - tol) && (s1_max < s2_max + tol));
}
  
  
// ----------------------------------------------------------------------------
inline vector<double> reduce_to_1D(const vector<Point2D>& points2D,
                                   Point2D& origin,
                                   Point2D& dir,
                                   const double tol)
// ----------------------------------------------------------------------------
{
  // Points are assumed to be collinear. Take first point as origin, and
  // distance to this point as the position along the axis (sign determined by
  // comparing with the next point that is not coincidental with the first)
  assert(points2D.size() > 1);
  vector<double> result(points2D.size(), 0);

  // Taking origin to be the first point
  origin = points2D[0];

  // Determining the direction by finding the next point that is not
  // coincidental with the origin
  uint ix = 1;
  while ( (dist(points2D[ix], origin) < tol) && ix != points2D.size())
    ++ix;
  if (ix == points2D.size()) {
    // all points coincidental to origin; we do not need to establish a
    // direction
    dir = {1, 0}; // just to avoid having it unitialized
    return result;
  }
    
  dir = points2D[ix] - origin;
  dir /= norm(dir);

  transform(points2D.begin(), points2D.end(), result.begin(),
                 [&origin, &dir] (const Point2D& p) {
                   const double val = sqrt(dist2(p, origin));
                   const Point2D d = p - origin;
                   return (d * dir > 0) ? val : -val;
                 });
  
  return result;
}

// ----------------------------------------------------------------------------
inline vector<Point2D>
reduce_to_2D(const vector<Point3D>& points3D,
             Point3D& origin, Point3D& v1, Point3D& v2, const double tol)
// ----------------------------------------------------------------------------
{
  // Points are assumed to be complanar. Take first point as origin, normalized
  // vector from first to second point to be the first axis, and normalized
  // vector, and the second axis to be obtained from the (negative) cross
  // product of the first axis and the plane normal.
  assert(points3D.size() > 2);

  // origin is the first point
  origin = points3D[0];
 
  // search for first axis
  uint ix = 1;
  Point3D n;
  while( (dist(origin, points3D[ix]) < tol) && ix < points3D.size())
    ++ix;
  if (ix == points3D.size()) {
    // points are all coincidental to origin.  We can choose any axes we want
    v1 = {1, 0, 0};
    v2 = {0, 1, 0};
    n  = {0, 0, 1};
  } else {
    v1 = (points3D[ix] - origin);
    v1 /= norm(v1);
    // determining second axis
    ++ix;
    Point3D tmp;
    while (ix != points3D.size() && (norm((points3D[ix] - origin) ^ v1) < tol))
      ++ix;
    if (ix == points3D.size()) {
      // points are collinear.  We can choose the second axis as we want, as
      // long as it does not coincide with the first.
      tmp = {v1[1], v1[2], v1[0]}; // rotate coordinates to make sure it is different
    } else {
      tmp = (points3D[ix] - origin);
    }
    n = v1 ^ tmp;
    n /= norm(n);
    v2 = n ^ v1; // computing v1
  }
  
  // we now have the orthonormal system (v1, v2, n)
  vector<Point2D> result(points3D.size());
  transform(points3D.begin(), points3D.end(), result.begin(),
                 [v1, v2, n, origin] (Point3D p) {
                   p -= origin;
                   Point3D tmp = solve_3D_matrix(v1, v2, n, p);
                   return Point2D {tmp[0], tmp[1]};
                 });
  return result;
}
  
// ----------------------------------------------------------------------------
bool point_inside_triangle(const Point2D& p1, const Point2D* const tri, const double tol)
// ----------------------------------------------------------------------------      
{
  return (projected_distance_to_line_2D(p1, tri[0], tri[1]) < tol) &&
         (projected_distance_to_line_2D(p1, tri[1], tri[2]) < tol) &&
         (projected_distance_to_line_2D(p1, tri[2], tri[0]) < tol);
}

// ----------------------------------------------------------------------------
inline bool bounding_boxes_overlap(const array<double, 4>& box1,
                                   const array<double, 4>& box2,
                                   const double tol)
// ----------------------------------------------------------------------------
{
  double dummy[2];
  return (isect_seg_seg_1D(&box1[0], &box2[0], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[2], &box2[2], tol, dummy) != DISJOINT);
}

// ----------------------------------------------------------------------------
inline bool bounding_boxes_overlap(const array<double, 6>& box1,
                                   const array<double, 6>& box2,
                                   const double tol)
// ----------------------------------------------------------------------------
{
  double dummy[2];
  return (isect_seg_seg_1D(&box1[0], &box2[0], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[2], &box2[2], tol, dummy) != DISJOINT) &&
         (isect_seg_seg_1D(&box1[4], &box2[4], tol, dummy) != DISJOINT);
}

// ----------------------------------------------------------------------------  
bool tri_inside_tri(const Point2D* tri1, const Point2D* tri2, const double tol)
// ----------------------------------------------------------------------------
{
  return point_inside_triangle(tri1[0], tri2, tol) &&
         point_inside_triangle(tri1[1], tri2, tol) &&
         point_inside_triangle(tri1[2], tri2, tol);
}

// ----------------------------------------------------------------------------  
array<Point2D, 2> scale_segment_beyond_bounding_box(const Point2D* seg,
                                              const array<double, 4>& bbox)
// ----------------------------------------------------------------------------
{
  const Point2D dp = seg[1] - seg[0];
  const uint comp = (fabs(dp[0]) > fabs(dp[1])) ? 0 : 1; // which is the largest component

  const double min_target = bbox[2*comp];
  const double max_target = bbox[2*comp + 1];

  array<Point2D, 2> result {seg[0], seg[1]};

  int expo = 1;
  
  while ( (min(result[0][comp], result[1][comp]) >= min_target) ||
          (max(result[0][comp], result[1][comp]) <= max_target)) {

    result[0] = seg[0] - (dp * pow(2, expo));
    result[1] = seg[1] + (dp * pow(2, expo));

    expo++;
  }
  return result;
}

// ----------------------------------------------------------------------------
enum ParallelState { IDENTICAL, PARALLEL, REGULAR };
ParallelState check_coplanar_state(const Point3D* const tri1,
                                   const Point3D* const tri2,
                                   double tol)
// ----------------------------------------------------------------------------
{
  const double ANGLE_TOL = 1e-6; // This is not a distance tolerance, so we don't use 'tol'
  Point3D n1; cross(tri1[1] - tri1[0], tri1[2] - tri1[0], n1);
  n1 /= norm(n1);
  Point3D n2; cross(tri2[1] - tri2[0], tri2[2] - tri2[0], n2);
  n2 /= norm(n2);

  Point3D tmp; cross(n1, n2, tmp);
  if (norm(tmp) < ANGLE_TOL) {
    // plane are parallel - but are they identical?  take the distance from one
    // point in each plane, and see if it has a non-neglible component along the
    // normal vector.
    Point3D dp = tri1[0] - tri2[0];
    const double sprod = dp * n1;
    return (fabs(sprod) < tol) ? IDENTICAL : PARALLEL;
  } 
  return REGULAR;
}


// ----------------------------------------------------------------------------
ParallelState check_colinear_state(const Point2D* const seg1,
                                   const Point2D* const seg2,
                                   double tol)
// ----------------------------------------------------------------------------
{
  const double ANGLE_TOL = 1e-6; // This is not a distance tolerance, so we don't use 'tol'
  Point2D u = seg1[1] - seg1[0]; u /= norm(u);
  Point2D v = seg2[1] - seg2[0]; v /= norm(v);
  const double tmp = u[0] * v[1] - u[1] * v[0];

  if (fabs(tmp) < ANGLE_TOL) {
    // lines defined by segments are parallel - but are they identical?
    Point2D w = seg2[0] - seg1[0];
    const double tmp2 = u[0] * w[1] - u[1] * w[0];

    return (fabs(tmp2) < tol) ? IDENTICAL : PARALLEL;
  }
  return REGULAR;
}

// ----------------------------------------------------------------------------
ParallelState check_colinear_state(const Point3D* const seg1,
                                   const Point3D* const seg2,
                                   double tol)
// ----------------------------------------------------------------------------
{
  const double ANGLE_TOL = 1e-6; // This is not a distance tolerance, so we don't use 'tol'
  Point3D u = seg1[1] - seg1[0]; u /= norm(u);
  Point3D v = seg2[1] - seg2[0]; v /= norm(v);
  const double tmp = norm(u ^ v);
  if (tmp < ANGLE_TOL) {
    // lines defined by segments are paralel - but are they colinear?
    Point3D w = seg2[0] - seg1[0];
    const double tmp2 = norm(u ^ w);
    return (tmp2 < tol) ? IDENTICAL : PARALLEL;
  }
  return REGULAR;
}

// ----------------------------------------------------------------------------
array<Point3D, 2> compute_plane_intersection(const Point3D* const tri1,
                                             const Point3D* const tri2)
// ----------------------------------------------------------------------------
{
  Point3D n1; cross(tri1[1] - tri1[0], tri1[2] - tri1[0], n1);
  n1 /= norm(n1);
  Point3D n2; cross(tri2[1] - tri2[0], tri2[2] - tri2[0], n2);
  n2 /= norm(n2);

  Point3D dir = n1 ^ n2; dir /= norm(dir); // direction vector

  // We now need to find a specific point on the intersection

  // determine index of largest component
  const auto m_el = max_element(&dir[0], &dir[0]+3,
                                [](double a, double b) {
                                  return fabs(a) < fabs(b);});
  const uint m_ix = uint(m_el - &dir[0]);
  const uint i1 = (m_ix + 1)%3;
  const uint i2 = (m_ix + 2)%3;

  const double d1 = tri1[0] * n1;
  const double d2 = tri2[0] * n2;
  
  const array<double, 2> rhs {d1, d2};

  // determining matrix column 1 for linear system to solve
  const array<double, 2> mcol1 {n1[i1], n2[i1]};
  const array<double, 2> mcol2 {n1[i2], n2[i2]};
  array<double, 2> res = solve_2D_matrix(mcol1, mcol2, rhs);
  
  Point3D p1;  p1[i1] = res[0]; p1[i2] = res[1]; p1[m_ix] = 0;
  return {p1, p1 + dir};
}

}; // end anonymous namespace


namespace TesselateUtils {

// ----------------------------------------------------------------------------
inline
IsectCase isect_point_seg_1D(const double p, const double* const s, const double tol)
// ----------------------------------------------------------------------------  
{
  const double smin = min(s[0], s[1]);
  const double smax = max(s[0], s[1]);

  if ((p < smin-tol) || (p > smax + tol)) return DISJOINT; // point outside bounding box
  if ((p < smin+tol) || (p > smax - tol)) return AT_POINT; // point "barely inside"

  return CONTAINING; // point well inside
  
}

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT or CONTAINING
IsectCase isect_point_seg_2D(const Point2D& p, const Point2D* const s, const double tol)
// ----------------------------------------------------------------------------  
{
  array<Point2D, 2> tmp {p, s[0]};
  if (check_colinear_state(&tmp[0], s, tol) == IDENTICAL) {
    Point2D origin, dir;
    const auto sys = reduce_to_1D({p, s[0], s[1]}, origin, dir, tol);

    return isect_point_seg_1D(sys[0], &sys[1], tol);
  }
  
  return DISJOINT;    
}
  
// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT or CONTAINING
IsectCase isect_point_seg_3D(const Point3D& p, const Point3D* const s, const double tol)
// ----------------------------------------------------------------------------  
{
  array<Point3D, 2> tmp{p, s[0]};
  if (check_colinear_state(&tmp[0], s, tol) == IDENTICAL) {
    Point3D origin, v1, v2;
    const auto sys = reduce_to_2D({p, s[0], s[1]}, origin, v1, v2, tol);
    
    return isect_point_seg_2D(sys[0], &sys[1], tol);
  }

  return DISJOINT;
}

// ----------------------------------------------------------------------------
inline  
IsectCase isect_seg_seg_1D(const double* const s1,
                            const double* const s2,
                            const double tol,
                            double* const result)
// ----------------------------------------------------------------------------  
{
  result[0] = nan(""); result[1] = nan(""); // value in case of no intersection
  const double s1_min = min(s1[0], s1[1]);
  const double s1_max = max(s1[0], s1[1]);
  const double s2_min = min(s2[0], s2[1]);
  const double s2_max = max(s2[0], s2[1]);

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
  result[0] = max(s1_min, s2_min);
  result[1] = min(s1_max, s2_max);
  return OVERLAPPING;
}
  
// ----------------------------------------------------------------------------
IsectCase isect_seg_seg_2D(const Point2D* const s1,
                           const Point2D* const s2,
                           const double tol,
                           Point2D* const result)
// ----------------------------------------------------------------------------    
{
  result[0] = {nan(""), nan("")};
  result[1] = {nan(""), nan("")};
  Point2D origin, dir; // will be used as output arguments further below
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_2D(s1, 2), bounding_box_2D(s2, 2), tol))
    return DISJOINT;

  // check if lines are colinear/parallel
  const auto clstate = check_colinear_state(s1, s2, tol);
  if (clstate == IDENTICAL) { 
    // segments lie on the same line.  Reduce dimensionality of problem
    const auto segs_1D = reduce_to_1D(vector<Point2D> {s1[0], s1[1], s2[0], s2[1]},
                                      origin, dir, tol);
    double res_1D[2];
    const auto icase = isect_seg_seg_1D(&segs_1D[0], &segs_1D[2], tol, res_1D);
    
    result[0] = origin + dir * res_1D[0];
    result[1] = origin + dir * res_1D[1];
    return (icase == OVERLAPPING) ? COLINEAR_OVERLAP : icase;
  } else if (clstate == PARALLEL) {
    // lines are parallel but not identical.  Cannot intersect
    return DISJOINT;
  }

  // Segments are not collinear. Compute unique intersection point of the two
  // lines on which segments lie.
  // We determine u and v such that: (1-u) * s1[0] + u * s1[1] = (1-v) *
  // s2[0] + v * s2[1]
  const auto uv = solve_2D_matrix( s1[1] - s1[0], s2[0] - s2[1], s2[0] - s1[0]);
  const double u = uv[0];
  
  result[0] = s1[0] * (1-u) + s1[1] * u; // this is the intersection point

  const auto pts_1D_1 = reduce_to_1D( vector<Point2D> {result[0], s1[0], s1[1]},
                                      origin, dir, tol);
  const auto icase1 = isect_point_seg_1D(pts_1D_1[0], &pts_1D_1[1], tol);

  const auto pts_1D_2 = reduce_to_1D( vector<Point2D> {result[0], s2[0], s2[1]},
                                      origin, dir, tol);
  const auto icase2 = isect_point_seg_1D(pts_1D_2[0], &pts_1D_2[1], tol);
  
  if ( (icase1 == DISJOINT) || (icase2 == DISJOINT))
    return DISJOINT; // intersection point clearly outside at least one segment

  if ( (icase1 == AT_POINT) || (icase2 == AT_POINT))
    return AT_POINT; // at least one of the segments intersect in an enpoint

  return OVERLAPPING; // if we got here, there was one intersection point that
                      // was clearly inside both segments
}

// ----------------------------------------------------------------------------
// Returns: DISJOINT, AT_POINT, COLINEAR_OVERLAP, CONTAINING or OVERLAPPING  
IsectCase isect_seg_seg_3D(const Point3D* const s1,
                           const Point3D* const s2,
                           const double tol,
                           Point3D* const result)
// ----------------------------------------------------------------------------    
{
  if (check_colinear_state(s1, s2, tol) == IDENTICAL) {
    Point3D origin, v1, v2;
    const auto sys = reduce_to_2D({s1[0], s1[1], s2[0], s2[1]}, origin, v1, v2, tol);
    array<Point2D, 2> res2D;
    const auto icase = isect_seg_seg_2D(&sys[0], &sys[2], tol, &res2D[0]);
    result[0] = origin + res2D[0][0] * v1 + res2D[0][1] * v2;
    result[1] = origin + res2D[1][0] * v1 + res2D[1][1] * v2;
    return icase;
  }
  return DISJOINT;
}
  
// ----------------------------------------------------------------------------    
IsectCase isect_seg_triangle_2D(const Point2D* const s, // pointer to two points
                                const Point2D* const tri, // pointer to three points
                                const double tol,
                                Point2D* const result)
// ----------------------------------------------------------------------------  
{
  result[0] = {nan(""), nan("")};
  result[1] = {nan(""), nan("")};
  
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_2D(s, 2), bounding_box_2D(tri, 3), tol)) 
    return DISJOINT;

  // Check if segment is fully contained within triangle.  We use negative
  // tolerance because if point touches boundary, we will capture it further
  // below.  (Otherwise, a boundary overlap would be interpreted as 'CONTAINING').
  if (point_inside_triangle(s[0], tri, tol) && point_inside_triangle(s[1], tri, -tol)) {
    result[0] = s[0];
    result[1] = s[1];
    return CONTAINING;
  }

  // Segment is not fully contained.  If an intersection exist, there has to be
  // an edge intersection somewhere.
  IsectCase icase[3];
  array<array<Point2D, 2>, 3> tmp_res;
  for (uint i = 0; i != 3; ++i) {
    array<Point2D, 2> tseg {tri[i], tri[(i+1)%3]};
    icase[i] = isect_seg_seg_2D(s, &tseg[0], tol, &(tmp_res[i][0]));

    if (icase[i] == CONTAINING || icase[i] == COLINEAR_OVERLAP) {
      result[0] = tmp_res[i][0];
      result[1] = tmp_res[i][1];
      return ALONG_BOUNDARY;
    }
  }
  // if we got here, there as no colinear overlap between segments.  There may
  // still be AT_POINT or regular OVERLAPs between segments

  vector<uint> overlaps;
  vector<uint> atpoints;
  for (uint i = 0; i!= 3; ++i) 
    if (icase[i] == OVERLAPPING) overlaps.push_back(i);
    else if (icase[i] == AT_POINT) atpoints.push_back(i);

  // we will return up to two intersection points, and if there are any overlap
  // points, we want to be sure that they are among the points returned, so we
  // add them first. (There may be two AT_POINTs that represent the same
  // triangle corner, in addition to an eventual overlap point, so if we do not
  // do this, we risk returning the two (identical) AT_POINTs).
  vector<Point2D> ipoints;
  for (uint i = 0; i != overlaps.size(); ++i) ipoints.push_back(tmp_res[overlaps[i]][0]);
  for (uint i = 0; i != atpoints.size(); ++i) ipoints.push_back(tmp_res[atpoints[i]][0]);

  for (uint i = 0; i != min((uint)ipoints.size(), (uint)2); ++i)
    result[i] = ipoints[i];

  return (overlaps.size() > 0) ? OVERLAPPING :
         (atpoints.size() > 0) ? AT_POINT :
                                 DISJOINT;
}

// ----------------------------------------------------------------------------
IsectCase isect_line_triangle_2D(const Point2D* const s, // segment defining the line
                                 const Point2D* const tri,
                                 const double tol,
                                 Point2D* const result)
// ----------------------------------------------------------------------------  
{
  // determining a segment large enough that it plays the role of an infinite
  // line in the segment/triangle intersection routine
  const auto longseg = scale_segment_beyond_bounding_box(s, bounding_box_2D(tri, 3));
  return isect_seg_triangle_2D(&longseg[0], tri, tol, result);
}

// ----------------------------------------------------------------------------  
IsectCase isect_triangle_triangle_2D(const Point2D* const tri1,
                                      const Point2D* const tri2,
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
  array<Point2D, 2> res_tmp; // we throw away the result for now
  for (uint i = 0; i != 3; ++i) {
    array<Point2D, 2> seg {tri1[i], tri1[(i+1)%3]};
    icase[i] = isect_seg_triangle_2D(&seg[0], tri2, tol, &res_tmp[0]);
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
IsectCase isect_triangle_triangle_3D(const Point3D* const tri1,
                                     const Point3D* const tri2,
                                     const double tol)
// ----------------------------------------------------------------------------    
{
  // check if bounding boxes overlap
  if (!bounding_boxes_overlap(bounding_box_3D(tri1, 3), bounding_box_3D(tri2, 3), tol))
    return DISJOINT;

  // check if triangles are coplanar
  const auto coplanar_state = check_coplanar_state(tri1, tri2, tol);
  if (coplanar_state == IDENTICAL) {
    Point3D origin, v1, v2; // we will not use these, but function requires them
    const auto system_2D = reduce_to_2D({tri1[0], tri1[1], tri1[2],
                     tri2[0], tri2[1], tri2[2]}, origin, v1, v2, tol);
    
    return isect_triangle_triangle_2D(&system_2D[0], &system_2D[3], tol);
    
  } else if (coplanar_state == PARALLEL) {
    return DISJOINT; // parallel but separate planes cannot intesect
    
  }

  // if we got here, we are in the 'regular' case where triangle planes do intersect
  
  // compute intersection line
  array<Point3D, 2> iline = compute_plane_intersection(tri1, tri2);

  // Reduce system to 2D
  Point3D origin, v1, v2;
  array<Point2D, 2> res_2D;
  array<Point3D, 2> res1;
  array<Point3D, 2> res2;

  const auto sys1_2D = reduce_to_2D({iline[0], iline[1], tri1[0], tri1[1], tri1[2]},
                                    origin, v1, v2, tol);
  const auto icase1 = isect_line_triangle_2D(&sys1_2D[0], &sys1_2D[2], tol, &res_2D[0]);
  res1[0] = origin + v1 * res_2D[0][0] + v2 * res_2D[0][1];
  res1[1] = origin + v1 * res_2D[1][0] + v2 * res_2D[1][1];

  const auto sys2_2D = reduce_to_2D({iline[0], iline[1], tri2[0], tri2[1], tri2[2]},
                                    origin, v1, v2, tol);
  const auto icase2 = isect_line_triangle_2D(&sys2_2D[0], &sys2_2D[2], tol, &res_2D[0]);
  res2[0] = origin + v1 * res_2D[0][0] + v2 * res_2D[0][1];
  res2[1] = origin + v1 * res_2D[1][0] + v2 * res_2D[1][1];

  if (icase1 == DISJOINT || icase2 == DISJOINT) // at least one triangle doesnt intersect at all
    return DISJOINT;

  if (icase1 == AT_POINT && icase2 == AT_POINT) // both triangles intersect line in a point
    return ((dist(res1[0], res2[0]) < tol) ||
            (dist(res1[0], res2[1]) < tol) ||
            (dist(res1[1], res2[0]) < tol) ||
            (dist(res1[1], res2[1]) < tol)) ? AT_POINT : DISJOINT;

  if (icase1 == AT_POINT) // if we got here, we know that icase2 intersects along a segment
    return (isect_point_seg_3D(res1[0], &res2[0], tol) == DISJOINT) ? DISJOINT :
                                                                      AT_POINT;
  if (icase2 == AT_POINT) // if we got here, we know that icase1 intersects algong a segment
    return (isect_point_seg_3D(res2[0], &res1[0], tol) == DISJOINT) ? DISJOINT :
                                                                  AT_POINT;
  // if we got here, both triangles intersect along a segment
  array<Point3D, 2> res_tmp;
  const auto icase3 = isect_seg_seg_3D(&res1[0], &res2[0], tol, &res_tmp[0]);
  return (icase3 == DISJOINT)                              ? DISJOINT :
         (icase3 == AT_POINT)                              ? AT_POINT :
    (icase1 == ALONG_BOUNDARY && icase2 == ALONG_BOUNDARY) ? ALONG_BOUNDARY : OVERLAPPING;
  //(icase1 == ALONG_BOUNDARY || icase2 == ALONG_BOUNDARY) ? ALONG_BOUNDARY : OVERLAPPING;
}

  
}; // end namespace TesselateUtils
