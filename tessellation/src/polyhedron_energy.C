#include "polyhedral_energies.h"
#include "interpoint_distances.h"
#include "tesselate_utils.h"
#include "common_defs.h"

using namespace std;
using namespace TesselateUtils;

namespace {

// // first entry in result array is the energy, the second the derivative
array<double, 2> energy(double dist, double R);

ValAndDer<Point3D> internal_energy(const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist);
  
ValAndDer<Point3D> boundary_energy(const Point3D* const bpoints,
                                   const uint num_bpoints,
                                   const Triangle* const btris,
                                   const unsigned int num_btris,
                                   const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist,
                                   const ClippedGrid<3>* const cgrid);

void add_boundary_contribution(const Point3D* const bpoints,
			       const Triangle& btri,
			       const Point3D* const ipoints,
			       const unsigned int num_ipoints,
                               const ClippedDomainType* const point_status,
			       const double vdist,
                               ValAndDer<Point3D>& result);

void accumulate_energy(const double dist,
		       const Point3D& p1,
		       const Point3D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point3D& p1_der_acc,
		       Point3D& p2_der_acc,
                       const bool are_mirror_points = false);


void add_outside_penalty_energy(const uint ipoint_ix,
                                const Point3D* const ipoints,
                                const Point3D* const bpoints,
                                const uint num_bpoints,
                                const Triangle* const btris,
                                const uint num_btris,
                                const double vdist,
                                ValAndDer<Point3D>& result);

// ValAndDer<Point3D> boundary_energy_2(const Point3D* const bpoints,
//                                      const uint num_bpoints,
//                                      const Triangle* const btris,
//                                      const unsigned int num_btris,
//                                      const Point3D* const ipoints,
//                                      const unsigned int num_ipoints,
//                                      const double vdist);

ClippedDomainType point_domain_type(const Point3D& pt,
                                    const Point3D* const bpoints,
                                    const Triangle* btris,
                                    const uint num_btris,
                                    const ClippedGrid<3>& cgrid);
  
}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer<Point3D> polyhedron_energy(const Point3D* const bpoints,
                                     const unsigned int num_bpoints,
                                     const Triangle* const btris,
                                     const unsigned int num_btris,
                                     const Point3D* const ipoints,
                                     const unsigned int num_ipoints,
                                     const double vdist,
                                     const ClippedGrid<3>* const cgrid)
// ----------------------------------------------------------------------------
{
  
  //compute internal energy (potential energy between internal points)
  const ValAndDer<Point3D> E_int = internal_energy(ipoints, num_ipoints, vdist);

  // // compute boundary energy (energy from interaction between internal points
  // // and boundary points, and from internal points and their mirror points)

  // compute boundary energy.  Divide 'vdist' by two to allow points to get a
  // bit closer to boundary than to each other (and tighten penalty for getting
  // even closer or trespassing boundary).
  const ValAndDer<Point3D> E_bnd = boundary_energy(bpoints, num_bpoints, btris,
                                                   num_btris, ipoints, num_ipoints,
                                                   vdist/1.5, cgrid); 

  // Adding up components and returning results
  ValAndDer<Point3D> E_tot = E_int;

  const double BND_FAC = 2; // if number of internal points and 'vdist' is well
                        // adjusted, this penalty factor should prevent
                        // trespassing of outer boundaries.
  E_tot.val += BND_FAC * E_bnd.val; 
  for (uint i = 0; i != (uint)E_tot.der.size(); ++i)
    E_tot.der[i] += BND_FAC * E_bnd.der[i];

  // E_tot.val += E_bnd.val; 
  // E_tot.der += E_bnd.der;
  
  return E_tot;
}

}; // end namespace TesselateUtils


namespace {

// ----------------------------------------------------------------------------
ValAndDer<Point3D> internal_energy(const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(ipoints, num_ipoints, vdist, false);

  // accumulating energies and storing total value and partial derivatives in 'result'
  ValAndDer<Point3D> result {0, vector<Point3D>(num_ipoints, {0.0, 0.0, 0.0})};
  for (const auto& d : dists)
    accumulate_energy(d.dist, ipoints[d.p1_ix], ipoints[d.p2_ix], vdist, 
		      result.val, result.der[d.p1_ix], result.der[d.p2_ix]);

  return result;
}

// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
		       const Point3D& p1,
		       const Point3D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point3D& p1_der_acc,
		       Point3D& p2_der_acc,
                       const bool are_mirror_points)
// ----------------------------------------------------------------------------  
{
  const array<double,2> e = energy(dist, vdist);

  // accumulating total energy
  energy_acc += e[0];
  
  // adding contribution to partial derivatives for the two points involved.
  // @@ Creation of temporary object on the line below.  If bottleneck, should
  // be rewritten.

  //If the points are mirror points of each other, all contributions to partial
  //derivatives should be doubled (and the accumulated values for the mirror
  //points are irrelevant)
  const double fac = are_mirror_points ? 2 : 1;
  
  const Point3D dvec = (p2 - p1) * (fac * e[1] / dist);
  p1_der_acc -= dvec;
  p2_der_acc += dvec;
  
}

// // ----------------------------------------------------------------------------
// uint dist_to_closest_tri(const Point3D& pt,
//                          const Point3D* const bpoints,
//                          uint num_bpoints,
//                          const Triangle* const btris,
//                          const unsigned int num_btris,
//                          Point3D& dir,
//                          double& dist)
// // ----------------------------------------------------------------------------  
// {
//   // find closest bpoint
//   const auto it = min_element(bpoints, bpoints + num_bpoints,
//                               [&pt] (const Point3D& p1, const Point3D& p2) {
//                                 return (dist2(p1, pt) < dist2(p2, pt));
//                               });
//   const uint pt_ix = uint(it - bpoints);

//   // find all triangles sharing this point  @@ Large potential for optimization!
//   vector<uint> neigh_tris;
//   for (uint i = 0; i != num_btris; ++i) 
//     if (find(btris[i].begin(), btris[i].end(), pt_ix) != btris[i].end())
//       neigh_tris.push_back(i);

//   // compute distances to the plane of these triangles
//   assert(neigh_tris.size() > 0);
//   array<Point3D, 3> tricorners;
//   dist = numeric_limits<double>::infinity();
//   uint closest_tri_ix = -1;
//   for (uint i = 0; i != (uint)neigh_tris.size(); ++i) {
//     const auto n = btris[neigh_tris[i]];
//     for (uint i = 0; i != 3; ++i)
//       tricorners[i] = bpoints[n[i]];

//     Point3D cur_dir;
//     double cur_dist = projected_distance_to_plane(pt, &tricorners[0], dir);
//     if (cur_dist < dist) {
//       dist = cur_dist;
//       dir = cur_dir;
//       closest_tri_ix = i;
//     }
//   }
//   return closest_tri_ix;
// }
  
// // ----------------------------------------------------------------------------
// ValAndDer<Point3D> boundary_energy_2(const Point3D* const bpoints,
//                                      const uint num_bpoints,
//                                      const Triangle* const btris,
//                                      const unsigned int num_btris,
//                                      const Point3D* const ipoints,
//                                      const unsigned int num_ipoints,
//                                      const double vdist)
// // ----------------------------------------------------------------------------
// {
//   ValAndDer<Point3D> result {0, vector<Point3D>(num_ipoints, {0.0, 0.0, 0.0})};
//   for (uint i = 0; i != num_ipoints; ++i) {
//     double dist;
//     Point3D dir; // normalized direction vector from point to closest point on triangle plane
//     dist_to_closest_tri(ipoints[i], bpoints, num_bpoints, btris, num_btris, dir, dist);
    
//     if (dist < 0) { // we are inside.  Compute energy as normal

//       const array<double, 2> e = energy(fabs(dist), vdist);
//       result.val += e[0];
//       result.der[i]  += dir * (-e[1]);
      
//     } else { // we are outside

//       const array<double, 2> e = energy(0, vdist);
//       result.val += e[0]  - dist * e[1];
//       result.der[i] -= dir * e[1];
//     }
//   }
//   return result;
// }

// ----------------------------------------------------------------------------
ClippedDomainType point_domain_type(const Point3D& pt,
                                    const Point3D* const bpoints,
                                    const Triangle* btris,
                                    const uint num_btris,
                                    const ClippedGrid<3>& cgrid)
// ----------------------------------------------------------------------------
{
  array<uint, 3> ix;
  for (uint i = 0; i != 3; ++i)
    if ((pt[i] < cgrid.bbox[2*i]) || (pt[i] > cgrid.bbox[2*i+1]))
      return OUTSIDE;
    else
      ix[i] = min((uint)max(floor((pt[i] - cgrid.bbox[2*i])/cgrid.cell_len[i]), 0.0),
                  cgrid.res[i] - 1);

  const auto type = cgrid.type[ix[0] + cgrid.res[0] * (ix[1] + cgrid.res[1] * ix[2])];

  if (type == FAR_INSIDE || type == CLOSE_INSIDE || type == OUTSIDE)
    return type;

  // if we got here, type is either UNDETERMINED or INTERSECTED. We must do an
  // explicit computation to determine its status vis-a-vis the boundary.
  bool on_bnd = false;
  const bool inside = inside_shell(pt, bpoints, btris, num_btris, 0, on_bnd);
  return ((!inside) && (!on_bnd)) ? OUTSIDE : UNDETERMINED;
}
  
// ----------------------------------------------------------------------------
ValAndDer<Point3D> boundary_energy(const Point3D* const bpoints,
                                   const uint num_bpoints,
                                   const Triangle* const btris,
                                   const unsigned int num_btris,
                                   const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist,
                                   const ClippedGrid<3>* const cgrid)
// ----------------------------------------------------------------------------
{
  ValAndDer<Point3D> result {0, vector<Point3D>(num_ipoints, {0.0, 0.0, 0.0})};

  // detect points lying outside
  vector<ClippedDomainType> point_status(num_ipoints);

  if (cgrid == nullptr) {
    // brute-force method that only distiguishes between points clearly outside
    // the shell and other points
    transform(ipoints, ipoints + num_ipoints, point_status.begin(), [&] (const Point3D& p) {
        bool on_bnd = false;
        bool inside = inside_shell(p, bpoints, btris, num_btris, 0, on_bnd);
        return ((!inside) && (!on_bnd)) ? OUTSIDE : UNDETERMINED;
    });
  } else {
    // benefit from the precomputed information found in 'cgrid' to
    // significantly speed up computations
    transform(ipoints,
              ipoints + num_ipoints,
              point_status.begin(),
              [&] (const Point3D& p) {
                return point_domain_type(p, bpoints, btris, num_btris, *cgrid);
              });
  }
    
  // looping across boundary faces and adding their energy contributions
  for (uint i = 0; i != num_btris; ++i) 
    add_boundary_contribution(bpoints, btris[i], ipoints, num_ipoints,
                              &point_status[0], vdist, result);

  
  // adding energy for points having exited the domain
  for (uint i = 0; i != num_ipoints; ++i) {
    if (point_status[i] == OUTSIDE)
      add_outside_penalty_energy(i, ipoints, bpoints, num_bpoints, btris,
                                 num_btris, vdist, result);
  }
  
  return result;
}

// ----------------------------------------------------------------------------  
void add_outside_penalty_energy(const uint ipoint_ix,
                                const Point3D* const ipoints,
                                const Point3D* const bpoints,
                                const uint num_bpoints,
                                const Triangle* const btris,
                                const uint num_btris,
                                const double vdist,
                                ValAndDer<Point3D>& result)
// ----------------------------------------------------------------------------
{
  // find the closest triangle the outside point projects to
  const double NUMTOL = sqrt(numeric_limits<double>::epsilon());
  const Point3D& p = ipoints[ipoint_ix];
  uint tri_ix;
  const Point3D cp = closest_point_on_triangle_surface(p, bpoints, num_bpoints,
                                                       btris, num_btris, tri_ix);

  // computing the derivative to be used
  const auto e = energy(0, vdist);
  const double der = fabs(e[1]); // derivative
  Point3D d = p - cp;
  const double dist = norm(d);
  if (dist < NUMTOL * vdist) {
    // we are basically on the triangle.  use triangle normal as direction vector
    array<Point3D, 3> tripts {bpoints[btris[tri_ix][0]],
                              bpoints[btris[tri_ix][1]],
                              bpoints[btris[tri_ix][2]]};
    cross(tripts[1] - tripts[0], tripts[2] - tripts[0], d);
    d /= norm(d);
  } else {
    d /= dist; // normalize direction vector
  }
  // const double penalty = e[0] + dist * der;

  // result.val += penalty;
  // result.der[ipoint_ix] += d * der;

  const double penalty = 10 * (dist * der);
  result.val += penalty;
  result.der[ipoint_ix] += 10 * der * d;
  
}
  
// ----------------------------------------------------------------------------
void add_boundary_contribution(const Point3D* const bpts, // boundary points
			       const Triangle& btri,
			       const Point3D* const ipoints,
			       const unsigned int num_ipoints,
                               const ClippedDomainType* const point_status,
			       const double vdist,
                               ValAndDer<Point3D>& result)
// ----------------------------------------------------------------------------
{
  // identify the points within reach of the face
  const array<Point3D, 3> tripts { bpts[btri[0]], bpts[btri[1]], bpts[btri[2]]};
  const double NUMTOL = sqrt(numeric_limits<double>::epsilon());

  for (uint i = 0; i != num_ipoints; ++i) {
    if (point_status[i] == FAR_INSIDE)
      continue;
    
    double dist_to_plane;
    Point3D proj_pt, edge_dir;
    ProjectionType ptype;
    if (point_on_triangle(ipoints[i], &tripts[0], vdist, proj_pt,
                          dist_to_plane, ptype, edge_dir)) {
      // this point is sufficiently close to the triangle that there will be an
      // energy involved.
      Point3D dir = proj_pt - ipoints[i];
      const double dist = norm(dir);
      if (dist < NUMTOL * vdist) {
        // direction 'dir' is degenerate.  We are basically on the triangle.  Use
        // triangle normal as direction vector.
        cross(tripts[1]-tripts[0], tripts[2] - tripts[0], dir);
        dir /= norm(dir); // normal points out of plane, which is the same
                          // direction as from point to plane, assuming the
                          // point is on the 'inside'
      } else {
        dir /= dist; // normalize direction vector
      }
      const array<double, 2> e = energy(dist, vdist);
      result.val += e[0];
      result.der[i] -= (dir * e[1]);
      
    }
  }
  
}



// ----------------------------------------------------------------------------
void add_boundary_contribution_old(const Point3D* const bpts, // boundary points
			       const Triangle& btri,
			       const Point3D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
                               ValAndDer<Point3D>& result)
{
  // identify the points within reach of the face
  const array<Point3D, 3> tripts { bpts[btri[0]], bpts[btri[1]], bpts[btri[2]]};
  const auto npoints = extract_from_range(ipoints, num_ipoints,
                                          [&tripts, &vdist] (const Point3D& p) {
                                            return point_on_triangle(p,
                                                                     tripts[0],
                                                                     tripts[1],
                                                                     tripts[2],
                                                                     vdist);});
  const auto& neigh_ixs = npoints.second;
  const auto& neigh_pts = npoints.first;

  // setting up a result structure only for the neigh points
  ValAndDer<Point3D> result_local {0, vector<Point3D>(neigh_pts.size(),
                                                      {0.0, 0.0, 0.0})};

  // computing the energy (and associated partial derivatives) between the
  // neighbor points and the two boundary point
  Point3D dummy; // used to store values we do not need
  for (uint i = 0; i != 3; ++i) {// loop over triangle corners
    const auto dvec = interpoint_distances(&neigh_pts[0], (uint)neigh_pts.size(),
                                           tripts[i], vdist);
    for (const DistanceEntry& d : dvec)
      accumulate_energy(d.dist, neigh_pts[d.p1_ix], tripts[i], vdist,
                        result_local.val, result_local.der[d.p1_ix], dummy);
  }
  
  // remapping results from result_local to result.  Multiplying by 1/6, since
  // the contribution from each corner will be counted as many times as the
  // valence of that corner. We do not keep track of valences; 6 is the typical
  // value for a regular mesh, which is why we divide it here.  @@ If this
  // simplification leads to bad results close to high valence (or low valence)
  // points, consider a more rigorous treatment.

  result.val += result_local.val * 0.1667; // 1/6
  for (uint i = 0; i != neigh_ixs.size(); ++i)
    result.der[neigh_ixs[i]] += result_local.der[i] * 0.1667;
  
  // Now, we compute the energy between neighbor points and their mirror images.
  // First, we remove the points that do not project perpendicularly to the
  // surface of the triangle itself.
  double dist_dummy;
  int sign_dummy;
  const double TOL = 1e-5 * vdist; // @@ should be enough?
  const auto per = extract_from_range(&neigh_pts[0], (uint)neigh_pts.size(),
                                      [&] (const Point3D& p) {
                                        return projects_to_triangle(p,
                                                                    &tripts[0],
                                                                    TOL,
                                                                    dist_dummy,
                                                                    sign_dummy);});
  const auto& per_pts_ixs = per.second;
  const auto& per_pts     = per.first;
  result_local.reset((uint)per_pts.size()); // reinitialize local result structure
 
  // computing mirror points and their energy contributions
  const auto mpoints = mirror_points_3D(&per_pts[0], (uint)per_pts.size(), &tripts[0]);
  const auto dvec = interpoint_distances(&per_pts[0], (uint)per_pts.size(),
                                         &mpoints[0], (uint)mpoints.size(), vdist);
  for (const DistanceEntry& d : dvec)  
    accumulate_energy(d.dist, per_pts[d.p1_ix], mpoints[d.p2_ix], vdist,
                      result_local.val, result_local.der[d.p1_ix], dummy, true);
  
  // Remapping results from 'results_local' to 'result'
  result.val += result_local.val;
  for (uint i = 0; i != per_pts.size(); ++i) 
    result.der[neigh_ixs[per_pts_ixs[i]]] += result_local.der[i];
}
  
// ----------------------------------------------------------------------------
// Energy as a function of distance.  Its support remains within the support of R
array<double, 2> energy(double dist, double R)
// ----------------------------------------------------------------------------
{
  const double tmp = max(R-dist, double(0));
  const double tmp2 = tmp*tmp;
  const double R4 = R*R*R*R; // normalizing factor so that E(0) = 1
  return {tmp2 * tmp2/R4, -4 * tmp*tmp*tmp/R4}; // energy and derivative
}

}; // end namespace Go
