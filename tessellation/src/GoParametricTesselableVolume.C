#include <assert.h>
#include <algorithm>
#include "common_defs.h"
#include "distance_function.h"
#include "GoParametricTesselableVolume.h"
#include "tesselate_parametric_volume.h"
#include "tesselate_polyhedron.h"

using namespace TesselateUtils;
using namespace std;
using namespace Go;

namespace {

typedef GoParametricTesselableVolume::PointType PointType;
typedef GoParametricTesselableVolume::FaceType FaceType;
  
// ----------------------------------------------------------------------------  
vector<Point2D> compute_surface_parameters(const vector<PointType>& boundary,
                                           const shared_ptr<const ParamSurface> surf);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
vector<Point3D> compute_volume_parameters(const vector<Go::Point> pts,
                                          const shared_ptr<const ParamVolume> pvol);
// ----------------------------------------------------------------------------  

  
// ----------------------------------------------------------------------------
uint find_vertex_index(const vector<shared_ptr<Vertex>>& vertices,
                       const shared_ptr<Vertex> vertex);
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------  
vector<pair<uint, bool>>
construct_edge_index_vector(const vector<shared_ptr<ftEdge>>& edges,
                            const shared_ptr<Loop> loop);
// ----------------------------------------------------------------------------  
  
  
}; // end anonymous namespace

namespace TesselateUtils
{

// ----------------------------------------------------------------------------
// Additional constructor
template<> template<>
GoParametricTesselableVolume::TesselableVolume(ftVolume& fvol)
// ----------------------------------------------------------------------------
{
  const auto shell = fvol.getOuterShell();
  
  // Setting the parametric volume
  volume_ = fvol.getVolume();

  // setting the corners (vertices)
  const auto vertices = fvol.vertices(); // shared_ptrs to vertices
  transform(vertices.begin(), vertices.end(), back_inserter(corners_),
            [] (const shared_ptr<Vertex> v) { return v->getVertexPoint();});
  
  // setting the edges 
  auto unique_edges = shell->getUniqueInnerEdges();
  transform(unique_edges.begin(), unique_edges.end(), back_inserter(edges_),
            [&vertices] (const shared_ptr<ftEdge> e) {
              // @@ for the moment, we require that the ParamCurve is not clipped
              cout << "tMin: " << e->tMin() << " " << "startparam: " << e->geomCurve()->startparam() << endl;
              cout << "tMax: " << e->tMax() << " " << "startparam: " << e->geomCurve()->endparam() << endl;

              // @@ DO SOMETHING IF THE PARAMETER INTERVALS ARE DIFFERENT

              //assert(e->tMin() == e->geomCurve()->startparam());
              //assert(e->tMax() == e->geomCurve()->endparam());
                
              return EdgeType { e->geomCurve(),
                                find_vertex_index(vertices, e->getVertex(true)), 
                                find_vertex_index(vertices, e->getVertex(false))};});
      
  // setting the faces.  We assume here that there are no inner loops
  auto ft_faces = shell->allFaces();
  transform(ft_faces.begin(), ft_faces.end(), back_inserter(faces_),
            [&unique_edges] (const shared_ptr<ftSurface> fs) {
              //assert(fs->nmbBoundaryLoops() == 1); // not yet support for inner loops
              return FaceType {
                fs->surface(),
                construct_edge_index_vector(unique_edges,
                                            fs->getBoundaryLoop(0))};});
  
}

// ----------------------------------------------------------------------------

template <> array<uint, 2>
GoParametricTesselableVolume::edgeCornerIndices(uint edge_ix) const
// ----------------------------------------------------------------------------  
{
  const auto e = edges_[edge_ix];
  return {get<1>(e), get<2>(e)};
}

// ----------------------------------------------------------------------------  
template <> vector<uint>
GoParametricTesselableVolume::faceBoundaryPointIndices(uint face_ix) const
// ----------------------------------------------------------------------------
{
  const auto edge_ixs = faces_[face_ix].ix;
  vector<uint> result;

  for (auto e : edge_ixs) {
    const auto edge = edges_[e.first];

    // internal points
    vector<uint> ixs(edge_ipoints_[e.first].size(), 0);
    iota(ixs.begin(), ixs.end(), edge_ipoints_start_ixs_[e.first]);

    
    if (e.second) { // oriented in the correct way

      result.push_back(get<1>(edge));

    } else { // reversely oriented

      result.push_back(get<2>(edge));
      reverse(ixs.begin(), ixs.end());
    }
    
    result.insert(result.end(), ixs.begin(), ixs.end());
  }

  return result;  
}

// ----------------------------------------------------------------------------
template<> 
void GoParametricTesselableVolume::
compute_tesselation(const array<PointType, 2>& boundary,
                    const EdgeType& edge,
                    const double vdist,
                    vector<PointType>& ipoints)
// ----------------------------------------------------------------------------
{
  // lowering slightly 'vdist' along the edge curves may help lower the
  // possibility of 'impossible' boundary triangles when tesselating the
  // interior volume later. @@ (not proven, just experienced)
  const auto param = tesselateParametricCurve(get<0>(edge), 0.9 * vdist);
  ipoints.resize(param.size() - 2);
  transform(param.begin()+1, param.end()-1, ipoints.begin(), [&edge]
            (const double d) { return PointType{get<0>(edge)->point(d), d}; });
}

// ----------------------------------------------------------------------------
vector<PointType> compute_3D_points(const shared_ptr<const ParamSurface> surf,
                                    const vector<Point2D>& parpoints)
// ----------------------------------------------------------------------------  
{
  vector<PointType> result;
  Go::Point gopoint;
  transform(parpoints.begin(), parpoints.end(), back_inserter(result),
            [&] (const Point2D& par) {
              surf->point(gopoint, par[0], par[1]);
              return PointType {gopoint, {par[0], par[1], nan("")}};
            });
  return result;
}
                       
// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<PointType>& boundary,
                    const FaceType& face,
                    const double vdist,
                    vector<PointType>& ipoints,
                    vector<Triangle>& triangles)
// ----------------------------------------------------------------------------  
{
  // determine surface parameter values for boundary points
  const vector<Point2D> par = compute_surface_parameters(boundary, face.surf);

  // computing tesselation
  const Mesh2D m2d = tesselateParametricSurface(face.surf, &par[0],
                                                (uint)par.size(), vdist);

  // compute internal 3D points
  ipoints = compute_3D_points(face.surf,
                              vector<Point2D> (m2d.points.begin() + boundary.size(),
                                               m2d.points.end()));

  // vector<Point3D> krull; // @@@
  // for (uint i = 0; i != boundary.size(); ++i) {
  //   const Go::Point pos = boundary[i].pos;
  //   krull.push_back(Point3D {pos[0], pos[1], pos[2]});
  // }
  // for (uint i = 0; i != ipoints.size(); ++i) {
  //   const Go::Point pos = ipoints[i].pos;
  //   krull.push_back(Point3D {pos[0], pos[1], pos[2]});
  // }

  
  triangles = m2d.tris;
  // fix orientation of triangles if necessary
}

  
// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<PointType>& bpoints,
                    const vector<Triangle>& btris,
                    const VolumeType& volume,
                    const double vdist,
                    vector<PointType>& ipoints,
                    vector<Tet>& tets)
// ----------------------------------------------------------------------------  
{
  // There is currently no parameter-based volume tesselation - the tesselation
  // is done directly in 3D space, and the corresponding parameters computed in
  // a second step.  Unlike the curve and surface case, it is possible to
  // tesselate the entity directly in 3D space, as we are here tesselating a
  // manifold that is not embedded in a higher-dimensional space.  On the other
  // hand, it might be more computationally efficient to tesselate in parameter
  // space (as in the curve and surface case), since we would then not need to
  // compute the parameters post-hoc.  So this could be considered later.

  // First, convert the boundary points to Point3D so that the polyhedron
  // tesselation routine can be called.
  vector<Point3D> bp3D(bpoints.size());
  transform(bpoints.begin(), bpoints.end(), bp3D.begin(),
            [&volume] (const PointType& p) {
              return Point3D {p.pos[0], p.pos[1], p.pos[2]};});

  const Mesh3D m3D = tesselatePolyhedron3D(&bp3D[0],
                                           (uint)bp3D.size(),
                                           &btris[0],
                                           (uint)btris.size(),
                                           vdist);

  const uint N = (uint)m3D.points.size(); // number of interior points

  // converting interior points to Go::Point
  vector<Go::Point> ipoints_go(N);
  transform(m3D.points.begin(), m3D.points.end(), ipoints_go.begin(),
            [](const Point3D& p) {return Go::Point(p[0], p[1], p[2]);});

  // computing interior point parameters
  const vector<Point3D> par = compute_volume_parameters(ipoints_go, volume); 
  //const vector<Point3D> par(ipoints_go.size(), {0, 0, 0});
  
  ipoints.resize(N);
  for (uint i = 0; i != N; ++i) 
    ipoints[i] = PointType(ipoints_go[i], par[i]);

  tets = m3D.tets;

  // ipoints = bpoints;
}

  
};


namespace {

  
// ----------------------------------------------------------------------------
uint find_vertex_index(const vector<shared_ptr<Vertex>>& vertices,
                       const shared_ptr<Vertex> vertex)
// ----------------------------------------------------------------------------
{
  const auto ptr = find(vertices.begin(), vertices.end(), vertex);
  assert(ptr != vertices.end()); // it should be in there somewhere
  return uint(ptr - vertices.begin());
}

// ----------------------------------------------------------------------------  
vector<pair<uint, bool>>
construct_edge_index_vector(const vector<shared_ptr<ftEdge>>& edges,
                            const shared_ptr<Loop> loop)
// ----------------------------------------------------------------------------  
{
  vector<pair<uint, bool>> result;
  const auto loop_edges = loop->getEdges(); // vector of shared_ptr<ftEdgeBase>
  for (auto e : loop_edges) {
    pair<uint, bool> cur_elem;

    // find the current edge in the edge list
    bool found = false;
    bool found_was_twin = false;
    for (uint i = 0; (i != edges.size()) && !found; ++i) {
      if (edges[i] == e) {
        cur_elem.first = i;
        found = true;
        found_was_twin = false;
      } else if (edges[i]->twin() == e.get()) {
        cur_elem.first = i;
        found = true;
        found_was_twin = true;
      }
    }
    assert(found);

    // determine orientation
    const auto ft_e = e->geomEdge();
    bool is_reversed = found_was_twin ^ ft_e->isReversed();
    cur_elem.second = !is_reversed;

    result.push_back(cur_elem);
  }
    
  return result;
}

// ----------------------------------------------------------------------------  
vector<Point2D> compute_surface_parameters(const vector<PointType>& boundary,
                                           const shared_ptr<const ParamSurface> surf)
// ----------------------------------------------------------------------------
{
  vector<Point2D> result(boundary.size());
  const double surf_extent = surf->containingDomain().diagLength();
  const double EPS = sqrt(numeric_limits<double>::epsilon()) * surf_extent;
  const double TOL = surf_extent * 1e-2;
  
  transform(boundary.begin(),
            boundary.end(),
            result.begin(),
            [&] (const PointType& p) {
              Go::Point cpoint; // closest point (we do not use this)
              double u, v; // parameters for closest point
              double dist; // distance to closest point (we do not use this)
              surf->closestPoint(p.pos, u, v, cpoint, dist , EPS);
              cout << "dist: " << dist << endl;
              if (dist > TOL) {
                throw runtime_error("Could not exactly locate point on surface");
              }
              return Point2D {u, v};
            });
  return result;
}

// ----------------------------------------------------------------------------
vector<Point3D> compute_volume_parameters(const vector<Go::Point> pts,
                                          const shared_ptr<const ParamVolume> pvol)
// ----------------------------------------------------------------------------
{
  vector<Point3D> result(pts.size());

  const Array<double, 6> pspan = pvol->parameterSpan();
  const double extent = max({(pspan[1] - pspan[0]),
                             (pspan[3] - pspan[2]),
                             (pspan[5] - pspan[4])});
  const double EPS = extent * sqrt(numeric_limits<double>::epsilon());
  const double TOL = extent * 2e-3; // @@ too large?
  
  transform(pts.begin(), pts.end(), result.begin(),
            [&] (const Go::Point& p) {
              Go::Point cpoint; // closest point( we do not use this)
              double u, v, w; // parameters for closest point
              double dist; // distance to closest point (we do not use this)
              pvol->closestPoint(p, u, v, w, cpoint, dist, EPS);
              if (dist > TOL) {
                cout << "Dist: " << dist << " TOL: " << TOL << endl;
                throw runtime_error("Could not exactly locate point inside volume.");
              }
              return Point3D {u, v, w};
            });
  return result;
}

  
} // end anonymous namespace
