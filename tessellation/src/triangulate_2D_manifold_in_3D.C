#include "triangulate_domain.h"
#include "tesselate_utils.h"

#include "GoTools/parametrization/PrFastUnorganized_OP.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
// #include "GoTools/parametrization/PrPrmUniform.h"
// #include "GoTools/parametrization/PrPrmLeastSquare.h"
// #include "GoTools/parametrization/PrPrmEDDHLS.h"
// #include "GoTools/parametrization/PrPrmMeanValue.h"

#include <memory>
#include <map>

using namespace Go;
using namespace std;
using namespace TesselateUtils;

namespace {

  struct InternalEdge {
    Segment s;
    uint tri_ix_1;
    uint tri_ix_2;
    bool processed;
  };
  
  vector<Point2D> compute_2D_paramerization(const Point3D* const points,
                                            const uint num_bpoints,
                                            const uint tot_num_points,
                                            const double vdist,
                                            const Point3D* const bpoint_normals); 

  vector<Point2D> compute_boundary_param(const Point3D* const bpoints,
                                         const Point3D* const bpoint_normals,
                                         const uint num_bpoints);
                                         
  vector<double> compute_approx_angles(const Point3D* const bpoints,
                                       const Point3D* const bpoint_normals,
                                       const uint num_bpoints);
  vector<Point2D> compute_2D_segments(const Point3D* const bpoints,
                                      const uint num_bpoints,
                                      const vector<double>& angles);
  
  double compute_angle(const Point3D& p1, const Point3D& p2,
                       const Point3D& p3, const Point3D& normal);

  void optimize_triangulation(vector<Triangle>& tris, const Point3D* const points);
  vector<InternalEdge> map_internal_edges(const vector<Triangle>& tri);
  void flip_this_edge(vector<InternalEdge>& edges, uint e_ix, vector<Triangle>& tris);  
  bool should_be_flipped(const InternalEdge& edge,
                         const vector<Triangle>& tris,
                         const Point3D* const points);

  inline double pow2(double a) {return a*a;}
}; // end anonymous namespace


namespace TesselateUtils {

// ============================================================================  
std::vector<Triangle> triangulate_2D_manifold_in_3D(const Point3D* const points,
                                                    const uint num_bpoints,
                                                    const uint tot_num_points,
                                                    const double vdist,
                                                    const Point3D* const bpoint_normals,
                                                    const Point2D* const param)
// ============================================================================  
{
  // first, establish a reasonable 2D parameterization of the 3D points
  const auto uv = (param) ?
    vector<Point2D> {param, param + tot_num_points} : 
    compute_2D_paramerization(points, num_bpoints, tot_num_points,
                              vdist, bpoint_normals);
    
  // then triangulate them
  const uint num_ipoints = tot_num_points - num_bpoints;
  const auto bbox = bounding_box_2D(&uv[0], (uint)uv.size());
  const double L1 = bbox[1] - bbox[0];
  const double L2 = bbox[3] - bbox[2];
  const double Lmin = min(L1, L2);
  const double Lmax = max(L1, L2);
  const double new_vdist = 6 * 4 * Lmax / (max(num_ipoints, (uint)1) * sqrt(Lmin / Lmax));
      
  // const double new_vdist = num_ipoints > 0 ?
  //   (8 * 4 / sqrt(num_ipoints)) :
  //   1.1 * max(bbox[1] - bbox[0], bbox[3] - bbox[2]);

  auto tri = triangulate_domain(&uv[0], num_bpoints, tot_num_points, new_vdist);

  // Post-processing step to ensure triangulation is optimal also in 3D space.

  //@@ untested.  If non-delaunay triangulations result from the above
  //procedure, the below function call should fix the problem.  Uncomment an
  //test if the situation arises!

  optimize_triangulation(tri, points);

  return tri;
  
}  

}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
vector<Point2D> compute_2D_paramerization(const Point3D* const points,
                                          const uint num_bpoints,
                                          const uint tot_num_points,
                                          const double vdist,
                                          const Point3D* const bpoint_normals)
// ----------------------------------------------------------------------------
{
  // the PrFastUnorganized_OP object below assumes boundary points are listed
  // last, so we need to reverse our usual order of boundary points and interior
  // points here
  vector<Point3D> tmp_points(points + num_bpoints, points + tot_num_points);
  tmp_points.insert(tmp_points.end(), points, points + num_bpoints);

  const uint num_ipoints = tot_num_points - num_bpoints;
  shared_ptr<PrFastUnorganized_OP> 
    u_points(new PrFastUnorganized_OP(tot_num_points,
                                      num_ipoints,
                                      (double*) &tmp_points[0]));
  u_points->setRadius(vdist);
  u_points->useRadius();
  u_points->initNeighbours();
  
  // Parametrize the boundary
  Point3D normal = compute_polygon_normal(points, num_bpoints);
  uint nix = uint(min_element(&normal[0], &normal[0]+3) - &normal[0]);
  Point3D sel_axis = {0, 0, 0};
  sel_axis[nix] = 1;

  const vector<Point3D> tmp_bnd_pts(points, points + num_bpoints);
  const vector<Point2D> bnd_par = compute_boundary_param(&tmp_bnd_pts[0],
                                                         bpoint_normals,
                                                         (uint)tmp_bnd_pts.size());
  
  for (uint i = 0; i != num_bpoints; ++i) {
    u_points->setU(i + num_ipoints, bnd_par[i][0]);
    u_points->setV(i + num_ipoints, bnd_par[i][1]);
  }

  // PrParametrizeBdy pb;
  // pb.setParamKind(PrCHORDLENGTHBDY);
  // pb.attach(u_points);
  // pb.parametrize();

  //Parametrize the interior
  //PrPrmMeanValue pi;
  PrPrmShpPres pi;
  pi.attach(u_points);
  pi.parametrize();

  // returning result
  vector<Point2D> result(tot_num_points);
  for (uint i = 0; i != num_bpoints; ++i)
    // paramter values of boundary points
    result[i] = {u_points->getU(i+num_ipoints), u_points->getV(i+num_ipoints)};

  for (uint i = 0; i != num_ipoints; ++i)
    // parameter values for interior points
    result [i + num_bpoints] = {u_points->getU(i), u_points->getV(i)};
                 
  return result;  
}

// ----------------------------------------------------------------------------
vector<Point2D> compute_boundary_param(const Point3D* const bpoints,
                                       const Point3D* const bpoint_normals,
                                       const uint num_bpoints)
// ----------------------------------------------------------------------------
{
  // Compute approximate 2D parameterization of boundary points, that tries as
  // much as possible to reflect the angles and distances between the 3D
  // boundary points.
  
  // compute angles and ensure they sum up to 2 pi
  const vector<double> angles = compute_approx_angles(bpoints,
                                                      bpoint_normals, num_bpoints);

  // compute segment lengths and directions, and adjust to ensure endpoints meet
  const vector<Point2D> segs = compute_2D_segments(bpoints, num_bpoints, angles);

  // compute 2D coordinates by adding up the segments
  vector<Point2D> result(num_bpoints);
  result[0] = Point2D {0.0, 0.0};
  for (int i = 1; (uint)i != num_bpoints; ++i) 
    result[i] = result[i-1] + segs[i-1];

  return result;
  
}

// ----------------------------------------------------------------------------  
vector<double> compute_approx_angles(const Point3D* const bpoints,
                                     const Point3D* const bpoint_normals,
                                     const uint num_bpoints)
// ----------------------------------------------------------------------------    
{
  const double PI = 3.1415926536;    
  // const Point3D midpt = accumulate(pts3D.begin(),
  //                                  pts3D.end(),
  //                                  Point3D {0.0, 0.0, 0.0}) / (double) N;
  
  vector<double> angles(num_bpoints);
  const Point3D normal = compute_polygon_normal(bpoints, num_bpoints);
  for (int i = 0; i != (int)num_bpoints; ++i) {
    angles[i] = compute_angle(bpoints[(i-1+num_bpoints)%num_bpoints],
                              bpoints[i],
                              bpoints[(i+1)%num_bpoints],
                              bpoint_normals ? bpoint_normals[i] : normal);
  }
  const double initial_angle_sum = accumulate(angles.begin(), angles.end(), 0.0);
  transform(angles.begin(), angles.end(), angles.begin(),
            [&](double x) {return x * 2 * PI / initial_angle_sum;});

  // the angles in 'angles' should now sum up to 2 PI
  return angles;
}

// ----------------------------------------------------------------------------  
double compute_angle(const Point3D& p1, const Point3D& p2, const Point3D& p3,
                     const Point3D& normal)
// ----------------------------------------------------------------------------  
{
  Point3D v1 = p2 - p1;
  Point3D v2 = p3 - p2;

  // projecting vectors onto plane
  const double nnorm2 = norm2(normal);
  v1 = v1 - ((v1 * normal) * normal / nnorm2);
  v2 = v2 - ((v2 * normal) * normal / nnorm2);
  
  const double v1_norm = norm(v1);
  const double v2_norm = norm(v2);

  const Point3D n = v1 ^ v2;
  //const Point3D c =  r ^ v1; 

  const double sign = (n*normal > 0) ? 1 : -1;
  double arg = v1 * v2 / (v1_norm * v2_norm);
  arg = arg > 1 ? 1 :
        arg < -1 ? -1 :
        arg;
  const double result = sign * acos(arg);
  return result;
}

// ----------------------------------------------------------------------------
vector<Point2D> compute_2D_segments(const Point3D* const bpoints,
                                    const uint num_bpoints,
                                    const vector<double>& angles)
// ----------------------------------------------------------------------------
{
  const int N = num_bpoints;
  vector<Point2D> result(N);
  
  // First, compute segments that have the right orientation and length
  double cur_angle = 0;
  for (int i = 0; i != N; ++i) {
    const double len = sqrt(dist2(bpoints[i],
                                  bpoints[(i+1)%N]));
    cur_angle += angles[i];
    result[i] = {len * cos(cur_angle), len * sin(cur_angle)};
  }

  // Then, adjust segment lengths so that endpoints match in the end
  const double x_err = accumulate(result.begin(), result.end(), 0.0,
                                  [](double acc, const Point2D& p)
                                  {return acc + p[0];});
  const double y_err = accumulate(result.begin(), result.end(), 0.0,
                                  [](double acc, const Point2D& p)
                                  {return acc + p[1];});

  // distributing error evently across segments, to ensure a set of segments
  // that sum up to zero
  const double dx = x_err/N;
  const double dy = y_err/N;
  transform(result.begin(), result.end(), result.begin(),
            [dx, dy] (const Point2D& p)
            { return Point2D {p[0] - dx, p[1] - dy};});
  return result;
}

// ----------------------------------------------------------------------------
void optimize_triangulation(vector<Triangle>& tris, const Point3D* const points)
// Flip internal edges until triangulation can no longer be improved
// ----------------------------------------------------------------------------
{
  vector<InternalEdge> edges = map_internal_edges(tris);

  // flip non-optimal edges until there are no more non-optimal edges left
  while (true) {
    bool changes_made = false;
    for (uint i = 0; i != (uint)edges.size(); ++i) {
      auto& e = edges[i];
      if ((e.processed == false && should_be_flipped(e, tris, points))) {
        flip_this_edge(edges, i, tris);
        changes_made = true;
      } else {
        e.processed = true;
      }
    }
    if (!changes_made)
      break;
  }
}

// ----------------------------------------------------------------------------
inline uint nonedge_corner(const Segment& s, const Triangle& t)
// ----------------------------------------------------------------------------  
{
  // return the single corner of t that is not on the edge s (it is assumed that
  // s is an edge of t)
  if ((t[0] != s[0]) && (t[0] != s[1])) return t[0];
  if ((t[1] != s[0]) && (t[1] != s[1])) return t[1];
  return t[2];
}    
    


// ----------------------------------------------------------------------------  
void flip_this_edge(vector<InternalEdge>& edges, uint e_ix, vector<Triangle>& tris)
// ----------------------------------------------------------------------------
{
  auto& e = edges[e_ix];  // reference to edge in question
  uint remote_corners[2];
  uint remote_corners_ix[2];
  const Triangle t1 = tris[e.tri_ix_1]; // these are the "old" triangles, to be
  const Triangle t2 = tris[e.tri_ix_2]; // replaced (so we take copies)

  // find the two triangle corners that are not lying on the common edge
  remote_corners[0] = nonedge_corner(e.s, t1);
  remote_corners[1] = nonedge_corner(e.s, t2);
  remote_corners_ix[0] = uint(find(t1.begin(), t1.end(), remote_corners[0]) - t1.begin());
  remote_corners_ix[1] = uint(find(t2.begin(), t2.end(), remote_corners[1]) - t2.begin());
  
  // changing mesh topology
  tris[e.tri_ix_1] = Triangle { remote_corners[0],
                                t1[(remote_corners_ix[0] + 1)%3], remote_corners[1] };
  tris[e.tri_ix_2] = Triangle { remote_corners[1],
                                t2[(remote_corners_ix[1] + 1)%3], remote_corners[0] };

  // updating current edge
  e.s = {remote_corners[0], remote_corners[1]};
  if (e.s[0] > e.s[1]) swap(e.s[0], e.s[1]); // ensure correct order for segment
  e.processed = true; // just flipped, so should be OK

  // updating other concerned edges
  uint dummy[3]; 
  for (auto& other_e : edges) {
    uint* const tri_ix_to_change =
      ((other_e.tri_ix_1 == e.tri_ix_1) ||
       (other_e.tri_ix_1 == e.tri_ix_2))    ? &(other_e.tri_ix_1) :
      ((other_e.tri_ix_2 == e.tri_ix_1) ||
       (other_e.tri_ix_2 == e.tri_ix_2))    ? &(other_e.tri_ix_2) :
      nullptr;

    if (tri_ix_to_change) {
      Triangle tmp(tris[e.tri_ix_1]);
      sort(tmp.begin(), tmp.end());
      *tri_ix_to_change =
        (set_difference(tmp.begin(), tmp.end(),
                        other_e.s.begin(), other_e.s.end(), dummy) - dummy == 1) ?
        e.tri_ix_1 : e.tri_ix_2;
    }

    other_e.processed = (bool)tri_ix_to_change; 
    
  }
}

// ----------------------------------------------------------------------------
bool should_be_flipped(const InternalEdge& edge,
                       const vector<Triangle>& tris,
                       const Point3D* const points)
// ----------------------------------------------------------------------------
{
  // find the corner in the second triangle that is _not_ on the shared edge
  const Triangle& t1 = tris[edge.tri_ix_1];
  const Triangle& t2 = tris[edge.tri_ix_2];

  const uint t2_corner_ix = nonedge_corner(edge.s, t2);
  
  // check if the found corner of the second triangle is within the sphere
  // defined by the three corners of the first triangle
  const Point3D& P1 = points[t1[0]];
  const Point3D& P2 = points[t1[1]];
  const Point3D& P3 = points[t1[2]];

  Point2D mcol1 {0,0}, mcol2 {0, 0}, rhs {0, 0}; // matrix columns and right
                                                 // hand side of linear system
  for (int i = 0; i != 3; ++i) {
    mcol1[0] += (P2[i] - P1[i]) * (P1[i] - P2[i]);
    mcol1[1] += (P2[i] - P1[i]) * (P2[i] - P3[i]);
    mcol2[0] += (P3[i] - P1[i]) * (P1[i] - P2[i]);
    mcol2[1] += (P3[i] - P1[i]) * (P2[i] - P3[i]);
    rhs[0] += 0.5 * (pow2(P1[i]) - pow2(P2[i])) - P1[i] * (P1[i] - P2[i]);
    rhs[1] += 0.5 * (pow2(P2[i]) - pow2(P3[i])) - P1[i] * (P2[i] - P3[i]);
  }
  const Point2D uv = solve_2D_matrix(mcol1, mcol2, rhs);
  const Point3D center = P1 + uv[0] * (P2 - P1) + uv[1] * (P3 - P1);

  double R2 = dist2(center, P1); // squared radius of sphere

  return dist2(points[t2_corner_ix], center) < R2; // corner is within the
                                                   // sphere => edge is not
                                                   // optimal
}


// ----------------------------------------------------------------------------
vector<InternalEdge> map_internal_edges(const vector<Triangle>& tris)
// ----------------------------------------------------------------------------
{
  // map all edges to their respective triangles
  multimap<Segment, uint> m; // map segment to triangle indices
  for (uint t_ix = 0; t_ix != tris.size(); ++t_ix) {
    const Triangle t = tris[t_ix];
    for (uint i = 0; i!= 3; ++i) {
      const uint ix1 = t[i];
      const uint ix2 = t[(i+1)%3];
      const Segment s = (ix1 < ix2) ? Segment {ix1, ix2} : Segment {ix2, ix1};
      m.insert(pair<Segment, uint>(s, t_ix));
    }
  }
      
  // identify internal edges and return them
  vector<InternalEdge> result;
  for (auto m_it = m.begin(); m_it != m.end(); ++m_it) {
    const Segment& key = m_it->first;
    if (m.count(key) > 1) {
      // this is an internal edge
      assert(m.count(key) == 2);
      const uint t_ix_1 = m_it->second;
      const uint t_ix_2 = (++m_it)->second;
      result.push_back( {key, t_ix_1, t_ix_2, false} );
    }
  }
  return result;
}



}; // end anonymous namespace 
