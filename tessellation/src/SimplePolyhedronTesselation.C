#include <iostream> // for debugging
#include <stdexcept> // for debugging
#include <array>
#include "tesselate_polyhedron.h"
#include "SimplePolyhedronTesselation.h"
#include "fit_points_to_plane.h"

using namespace std;
using namespace TesselateUtils;

namespace {
  // return all the corners, in order
  vector<uint> assemble_segments(const vector<Segment>& segs);

  vector<Point3D> center_points(const vector<Point3D>& pts, Point3D& offset);
  array<Point3D, 2> compute_plane_coord_system(const Point3D& plane_normal);
  vector<Point2D> project_points_to_plane(const vector<Point3D>& pts,
                                          const array<Point3D, 2>& plane_axes);
  vector<Point3D> lift_plane_to_3D(const vector<Point2D>& pts, 
                                   const array<Point3D, 2>& plane_axes,
                                   const Point3D& offset);
}; // end anonymous namespace

namespace TesselateUtils
{

// ----------------------------------------------------------------------------  
template<>  
array<uint, 2> SimplePolyhedron::edgeCornerIndices(uint edge_ix) const
// ----------------------------------------------------------------------------  
{
  return edges_[edge_ix];
}

// ----------------------------------------------------------------------------  
template<>  
vector<uint> SimplePolyhedron::faceBoundaryPointIndices(uint face_ix) const
// ----------------------------------------------------------------------------
{
  // first, we identify the relevant corner points
  const auto face_edges = faces_[face_ix].edge_ixs;
  const vector<uint> corners = assemble_segments(extract_from_range(edges_, face_edges));

  vector<uint> result;
  uint corner_ix = 0;
  // then, we add the indices of the points internal to each edge
  for(const auto& e_ix : face_edges) {
    result.push_back(corners[corner_ix]);
    vector<uint> ixs(edge_ipoints_[e_ix].size(), 0);
    iota(ixs.begin(), ixs.end(), edge_ipoints_start_ixs_[e_ix]);

    // if edge is oriented the opposite direction, we need to flip the order of internal points
    if (edges_[e_ix][0] != corners[corner_ix]) {
      reverse(ixs.begin(), ixs.end());
    }
    
    result.insert(result.end(), ixs.begin(), ixs.end());
    ++corner_ix;
  }
  return result;
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const array<Point3D, 2>& boundary,
                                           const Segment& edge,
                                           const double vdist,
                                           vector<Point3D>& ipoints)
// ----------------------------------------------------------------------------
{
  vector<Point3D> points = tesselateSegment(boundary[0], boundary[1], vdist);
  ipoints.resize(points.size()-2);
  copy(points.begin() + 1, points.end()-1, ipoints.begin());
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const vector<Point3D>& boundary,
                                           const FaceLoop& face,
                                           const double vdist,
                                           vector<Point3D>& ipoints,
                                           vector<Triangle>& triangles)
// ----------------------------------------------------------------------------
{
  // fit a plane to the points on the boundary, with normal oriented such loop is
  // counterclockwise
  const array<double, 4> pl = fit_loop_to_plane(&boundary[0],
                                                (uint)boundary.size());

  // check if any point is too far off the fitted plane
  for (const auto& p : boundary) {
    const double dist = p[0] * pl[0] + p[1] * pl[1] + p[2] * pl[2] + pl[3];
    if (dist > vdist/2)
      throw runtime_error("The boundary points of the face are not sufficiently"
                          " close to a common plane.");
  }
  // if we got here, we know that the plane is OK
  Point3D offset;
  const vector<Point3D> centered_boundary = center_points(boundary, offset);
  const array<Point3D, 2> plane_axes = compute_plane_coord_system( {pl[0], pl[1], pl[2]} );
  const vector<Point2D> projpoints = project_points_to_plane(centered_boundary, plane_axes);

  // tesselate face in plane
  const Mesh2D m2d = tesselatePolygon2D(&projpoints[0],
                                        (uint) projpoints.size(),
                                        vdist, NULL);
  
  // lift interior points of face back up into 3D
  auto ipoints2D = vector<Point2D>(m2d.points.begin() + boundary.size(), m2d.points.end());
  ipoints = lift_plane_to_3D(ipoints2D, plane_axes, offset);
  triangles = m2d.tris;

  // fix orietation of triangles if necessary
  if (!face.ccw)
    transform(triangles.begin(), triangles.end(), triangles.begin(), [](const Triangle& t) {
        return Triangle {t[1], t[0], t[2]}; // flip orientation of triangle
      });
  
}

// ----------------------------------------------------------------------------  
template<> 
void SimplePolyhedron::compute_tesselation(const vector<Point3D>& bpoints,
                                           const vector<Triangle>& btris,
                                           const VolumeType& volume,  // dummy; empty type
                                           const double vdist,
                                           vector<Point3D>& ipoints,
                                           vector<Tet>& tets)
// ----------------------------------------------------------------------------
{
  const Mesh3D m3d = tesselatePolyhedron3D(&bpoints[0], (uint)bpoints.size(),
                                           &btris[0], (uint)btris.size(),
                                           vdist);
  ipoints = m3d.points;
  tets = m3d.tets;
}
  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
// return all the corners, in order
vector<uint> assemble_segments(const vector<Segment>& segs)
// ----------------------------------------------------------------------------
{
  vector<uint> result;
  if (!segs.empty()) {
    result.push_back(segs[0][0]);
    for (const auto s : segs)
      result.push_back(s[0] == result.back() ? s[1] : s[0]);
  }
  result.pop_back(); // remove duplicate point at end
  return result;
}

// ----------------------------------------------------------------------------  
vector<Point3D> center_points(const vector<Point3D>& pts, Point3D& offset)
// ----------------------------------------------------------------------------  
{
  offset = accumulate(pts.begin(), pts.end(), Point3D {0.0, 0.0}) / (double)pts.size();
  vector<Point3D> result(pts.size());
  transform(pts.begin(), pts.end(), result.begin(), [offset] (const Point3D& p) {
      return p-offset; });
  return result;
}

// ----------------------------------------------------------------------------  
array<Point3D, 2> compute_plane_coord_system(const Point3D& n)
// ----------------------------------------------------------------------------  
{
  const vector<Point3D> candidates {n ^ Point3D {1, 0, 0},
                                    n ^ Point3D {0, 1, 0},
                                    n ^ Point3D {0, 0, 1} };
  const Point3D v1 = *max_element(candidates.begin(),
                                  candidates.end(),
                                  [](const Point3D& p1, const Point3D& p2) {
                                    return norm2(p1) < norm2(p2);
                                  });
  const Point3D v2 = n ^ v1;

  return { v1 / sqrt(norm2(v1)), v2 / sqrt(norm2(v2)) };
  
}

// ----------------------------------------------------------------------------  
vector<Point2D> project_points_to_plane(const vector<Point3D>& pts,
                                        const array<Point3D, 2>& axes)
// ----------------------------------------------------------------------------  
{
  vector<Point2D> result(pts.size());
  transform(pts.begin(), pts.end(), result.begin(), [&axes] (const Point3D& p) {
      return Point2D { p * axes[0], p * axes[1]}; // components of p along the axes
    });
  return result;
}

// ----------------------------------------------------------------------------  
vector<Point3D> lift_plane_to_3D(const vector<Point2D>& pts, 
                                 const array<Point3D, 2>& axes,
                                 const Point3D& offset)
// ----------------------------------------------------------------------------  
{
  vector<Point3D> result(pts.size());
  transform(pts.begin(), pts.end(), result.begin(), [&axes, &offset] (const Point2D& p) {
      return (p[0] * axes[0]) + (p[1] * axes[1]) + offset;
    });
  return result;
}

  
}; // end anonymous namespace

