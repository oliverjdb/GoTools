#include "polyhedral_energies.h"
#include "interpoint_distances.h"
#include "tesselate_utils.h"
#include "common_defs.h"
#include <fstream> // @@ for debugging purposes

using namespace std;
using namespace TesselateUtils;

namespace {

// first entry in result array is the energy, the second the derivative
array<double, 2> energy(double dist, double R);

ValAndDer<Point2D> internal_energy(const Point2D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist);
                   
ValAndDer<Point2D> boundary_energy(const Point2D* const bpoints,
                                   const unsigned int num_bpoints,
                                   const Point2D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist,
                                   const ClippedGrid<2>* const cgrid);

void add_boundary_contribution(const Point2D& bp1,
                               const Point2D& bp2,
                               const Point2D* const ipoints,
                               const unsigned int num_ipoints,
                               const ClippedDomainType* const point_status, 
                               const double vdist,
                               ValAndDer<Point2D>& result);

// void add_boundary_contribution(const Point2D& bp1,
// 			       const Point2D& bp2,
// 			       const Point2D* const ipoints,
// 			       const unsigned int num_ipoints,
// 			       const double vdist,
//                                ValAndDer<Point2D>& result);

void accumulate_energy(const double dist,
		       const Point2D& p1,
		       const Point2D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point2D& p1_der_acc,
		       Point2D& p2_der_acc,
                       const bool are_mirror_points= false);
                       
void add_outside_penalty_energy(const uint ipoint_ix,
                                const Point2D* const ipoints,
                                const Point2D* const bpoints,
                                const uint num_bpoints,
                                const double vdist,
                                ValAndDer<Point2D>& result);

ClippedDomainType point_domain_type(const Point2D& pt,
                                    const Point2D* const bpoints,
                                    const uint num_bpoints,
                                    const ClippedGrid<2>& cgrid);
  
}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer<Point2D> polygon_energy(const Point2D* const bpoints,
                                  const unsigned int num_bpoints,
                                  const Point2D* const ipoints,
                                  const unsigned int num_ipoints,
                                  const double vdist,
                                  const ClippedGrid<2>* const cgrid)
// ----------------------------------------------------------------------------
{
  
  // compute internal energy (potential energy between internal points)
  const ValAndDer<Point2D> E_int = internal_energy(ipoints, num_ipoints, vdist);

  // compute boundary energy (energy from interaction between internal points
  // and boundary points, and from internal points and their mirror points)
  const ValAndDer<Point2D> E_bnd = boundary_energy(bpoints, num_bpoints,
                                                   ipoints, num_ipoints,
                                                   vdist/1.5, cgrid);

  // Adding up components and returning results
  ValAndDer<Point2D> E_tot = E_int;

  const double BND_FAC = 2; // 2;// Increase penalty for approaching border
  E_tot.val += BND_FAC * E_bnd.val;
  for (uint i = 0; i != (uint)E_tot.der.size(); ++i)
    E_tot.der[i] += E_bnd.der[i] * BND_FAC;
  
  //E_tot.der += E_bnd.der;
  
  return E_tot;
}

}; // end namespace Go


namespace {

// ----------------------------------------------------------------------------
ValAndDer<Point2D> internal_energy(const Point2D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(ipoints, num_ipoints, vdist, false);

  // accumulating energies and storing total value and partial derivatives in 'result'
  ValAndDer<Point2D> result {0, vector<Point2D>(num_ipoints, {0.0, 0.0})};
  for (const auto& d : dists)
    accumulate_energy(d.dist, ipoints[d.p1_ix], ipoints[d.p2_ix], vdist, 
		      result.val, result.der[d.p1_ix], result.der[d.p2_ix]);

  return result;
}

// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
		       const Point2D& p1,
		       const Point2D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point2D& p1_der_acc,
		       Point2D& p2_der_acc,
                       const bool are_mirror_points)
// ----------------------------------------------------------------------------  
{
  const array<double,2> e = energy(dist, vdist);

  // accumulating total energy
  energy_acc += e[0];
  
  // adding contribution to partial derivatives for the two points involved.  @@
  // Creation of temporary object on the line below.  If bottleneck, should be
  // rewritten.

  //If the points are mirror points of each other, all contributions to partial
  //derivatives should be doubled (and the accumulated values for the mirror
  //points are irrelevant)
  const double fac = are_mirror_points ? 2 : 1;

  // computing partial derivatives
  const Point2D dvec = (p2 - p1) * (fac * e[1] / dist);
  p1_der_acc -= dvec;
  p2_der_acc += dvec;
  
}
						 
// // ----------------------------------------------------------------------------
//   ValAndDer<Point2D> boundary_energy(const Point2D* const bpoints,
// 			  const unsigned int num_bpoints,
// 			  const Point2D* const ipoints,
// 			  const unsigned int num_ipoints,
// 			  const double vdist)
// // ----------------------------------------------------------------------------
// {
//     ValAndDer<Point2D> result {0, vector<Point2D>(num_ipoints, {0.0, 0.0})};

//   // looping across boundary segments and adding their energy contributions
//   for (uint i = 0; i != num_bpoints; ++i)
//     add_boundary_contribution(bpoints[i], bpoints[(i+1)%num_bpoints],
// 			      ipoints, num_ipoints, vdist, result);
//   return result;
// }

// ----------------------------------------------------------------------------  
inline ClippedDomainType point_domain_type(const Point2D& pt,
                                           const Point2D* const bpoints,
                                           const uint num_bpoints,
                                           const ClippedGrid<2>& cgrid)
// ----------------------------------------------------------------------------  
{
  const double& dx = cgrid.cell_len[0];
  const double& dy = cgrid.cell_len[1];
  const uint ix = min((uint)max(floor((pt[0] - cgrid.bbox[0])/dx), 0.0), cgrid.res[0]-1);
  const uint iy = min((uint)max(floor((pt[1] - cgrid.bbox[2])/dy), 0.0), cgrid.res[1]-1);

  const auto type = (cgrid.type[iy * cgrid.res[0] + ix]);

  if (type == FAR_INSIDE || type == CLOSE_INSIDE || type == OUTSIDE)
    return type;

  // if we got here, type is either UNDETERMINED or INTERSECTED.  We must do an
  // explicit computation to determine its status vis-a-vis the boundary.
  bool on_bnd = false;
  const bool inside = inpolygon(pt, bpoints, num_bpoints, 0, on_bnd);

  return (inside || on_bnd) ? CLOSE_INSIDE : OUTSIDE;
}
  
// ----------------------------------------------------------------------------
ValAndDer<Point2D> boundary_energy(const Point2D* const bpoints,
                                   const unsigned int num_bpoints,
                                   const Point2D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist,
                                   const ClippedGrid<2>* const cgrid)
// ----------------------------------------------------------------------------
{
  ValAndDer<Point2D> result {0, vector<Point2D>(num_ipoints, {0.0, 0.0})};
  
  // detect points lying outside
  vector<ClippedDomainType> point_status(num_ipoints);

  //if (true) {// (cgrid==nullptr) {
  if (cgrid==nullptr) {
    // simple, brute-force method that only distinguishes between points that
    // are clearly outside polygon and other points
    transform(ipoints, ipoints + num_ipoints, point_status.begin(), [&] (const Point2D& p) {
        bool on_bnd = false;
        bool inside = inpolygon(p, bpoints, num_bpoints, 0, on_bnd);

        return ((!inside) && (!on_bnd)) ? OUTSIDE : UNDETERMINED;
      });
  } else {
    // benefit from the precomputed information found in 'cgrid' to
    // significantly speed up computations
    transform(ipoints, ipoints + num_ipoints, point_status.begin(), [&] (const Point2D& p) {
        return point_domain_type(p, bpoints, num_bpoints, *cgrid);
      });
  }

  //@@
  // ofstream os_poly("poly.mat");
  // for (uint i = 0; i != num_bpoints; ++i) {
  //   os_poly << bpoints[i];
  // }
  // os_poly.close();
  // ofstream os_cgrid("cgrid.mat");
  // for (uint i = 0; i != cgrid->type.size(); ++i) {
  //   os_cgrid << cgrid->type[i] << " ";
  // }
  // os_cgrid.close();

  
  // looping across boundary segments and adding their energy contributions
  for (uint i = 0; i != num_bpoints; ++i)
    add_boundary_contribution(bpoints[i], bpoints[(i+1)%num_bpoints],
			      ipoints, num_ipoints, &point_status[0], vdist, result);

  // adding energy for points having exited the domain
  for (uint i = 0; i != num_ipoints; ++i) {
    if (point_status[i] == OUTSIDE)
      add_outside_penalty_energy(i, ipoints, bpoints, num_bpoints, vdist, result);
  }
  
  return result;
}
// ----------------------------------------------------------------------------
  void add_outside_penalty_energy(const uint ipoint_ix,
                                  const Point2D* const ipoints,
                                  const Point2D* const bpoints,
                                  const uint num_bpoints,
                                  const double vdist,
                                  ValAndDer<Point2D>& result)
// ----------------------------------------------------------------------------
{
  // find the closest point on the boundary, as well as the segment it lies on
  const double NUMTOL = sqrt(numeric_limits<double>::epsilon());
  const Point2D& p = ipoints[ipoint_ix];
  uint seg_ix;
  const Point2D cp = closest_point_on_loop(p, bpoints, num_bpoints, seg_ix);

  // computing the derivative to be used
  const auto e = energy(0, vdist);
  const double der = fabs(e[1]);
  Point2D d = p - cp;
  const double dist = norm(d);
  if (dist < NUMTOL * vdist) {
    // we are basically on the segment.  Use its outwards-pointing normal as direction vector
    const Point2D tangent = bpoints[(seg_ix+1)%num_bpoints] - bpoints[seg_ix];
    d = Point2D {tangent[1], -tangent[0]};
    d /= norm(d);
  } else {
    d /= dist;
  }
  const double penalty = e[0] + dist*der;
  result.val += penalty;
  result.der[ipoint_ix] += d * der;
  
}
  
// ----------------------------------------------------------------------------  
void add_boundary_contribution(const Point2D& bp1,
			       const Point2D& bp2,
			       const Point2D* const ipoints,
			       const unsigned int num_ipoints,
                               const ClippedDomainType* const point_status,
			       const double vdist,
			         ValAndDer<Point2D>& result)
// ----------------------------------------------------------------------------
{
  const double NUMTOL = sqrt(numeric_limits<double>::epsilon());
  for (uint i = 0; i != num_ipoints; ++i) {
    if ((point_status[i] == OUTSIDE) || (point_status[i] == FAR_INSIDE))
      continue;
    
    if (point_on_line_segment(ipoints[i], bp1, bp2, vdist, false)) {
      // if (point_status[i] == FAR_INSIDE) {
      //   throw runtime_error("this should not happen!");
      // }
      bool at_corner; // given a value in the function call below
      const Point2D proj_pt =
        projected_point_on_segment(ipoints[i], bp1, bp2, at_corner);
      // this point is sufficiently close to the line segment that there will be
      // an energy involved
      Point2D dir = proj_pt - ipoints[i];
      const double dist = norm(dir);
      if (dist < NUMTOL * vdist) {
        // direction 'dir' is degenerate.  We are basically _on_ the line
        // segment.  Use outward-pointing normal as direction vector
        Point2D tangent = bp2 - bp1;
        dir[0] = tangent[1];
        dir[1] = -tangent[0];
        dir /= norm(dir); 
      } else {
        dir /= dist;
      }
      const array<double, 2> e = energy(dist, vdist);
      result.val += e[0];
      result.der[i] -= (dir * e[1]);
    }
  }
}
  
// ----------------------------------------------------------------------------  
void add_boundary_contribution_old(const Point2D& bp1,
			       const Point2D& bp2,
			       const Point2D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
			         ValAndDer<Point2D>& result)
// ----------------------------------------------------------------------------
{
  // identify the points within reach of the boundary ("neighbor points"). 
  const auto npoints = extract_from_range(ipoints, num_ipoints,
  			     [&bp1, &bp2, vdist](const Point2D& p) {
        return point_on_line_segment(p, bp1, bp2, vdist, false);});
  const auto& neigh_ixs    = npoints.second;   
  const auto& neigh_points = npoints.first; 

  // Setting up a result structure only for the neigh points
    ValAndDer<Point2D> result_local {0, vector<Point2D>(neigh_points.size(), {0.0, 0.0})};
  Point2D dummy; // used to store values we don't want
  
  // compute the energy (and associated partial derivatives) between the
  // neighbor points and the two boundary points.  
  vector<DistanceEntry> dvec;
  dvec = interpoint_distances(&neigh_points[0], (uint)neigh_points.size(), bp1, vdist);
  for (const DistanceEntry& d : dvec)
    accumulate_energy(d.dist, neigh_points[d.p1_ix], bp1, vdist,
		      result_local.val, result_local.der[d.p1_ix], dummy);
  dvec = interpoint_distances(&neigh_points[0], (uint)neigh_points.size(), bp2, vdist);
  for (const auto& d : dvec)
    accumulate_energy(d.dist, neigh_points[d.p1_ix], bp2, vdist, 
		      result_local.val, result_local.der[d.p1_ix], dummy);

  // remapping results from result_local to result.  Multiplying by 0.5, since
  // segment endpoints will be counted again when neighbour segments are treated
  result.val += result_local.val * 0.5;
  for (uint i = 0; i != neigh_ixs.size(); ++i) 
    result.der[neigh_ixs[i]] += result_local.der[i] * 0.5;
  
  // Now, we compute the energy between the neighbor points and their mirror images.
  // first, remove points that do not project perpendicularly to within the segment.
  const auto per = extract_from_range(&neigh_points[0], (uint)neigh_points.size(),
				      [&bp1, &bp2] (const Point2D& p) {
					return projects_to_segment(p, bp1, bp2);});
  const auto& per_pts_ixs = per.second;
  const auto& per_pts     = per.first;
  result_local.reset((uint)per_pts.size()); // reinitialize local result structure

  const auto mpoints = mirror_points_2D(&per_pts[0], (uint)per_pts.size(), bp1, bp2);
  dvec = interpoint_distances(&per_pts[0], (uint)per_pts.size(),
			      &mpoints[0], (uint)per_pts.size(), vdist);

  for (const auto& d : dvec) 
      accumulate_energy(d.dist, per_pts[d.p1_ix], mpoints[d.p2_ix], vdist,
                        result_local.val, result_local.der[d.p1_ix], dummy, true);

  // remapping results from result_local to result
  result.val += result_local.val;
  for (uint i = 0; i != per_pts.size(); ++i)
    result.der[neigh_ixs[per_pts_ixs[i]]] += result_local.der[i];
}

  
// // ----------------------------------------------------------------------------
// // Energy as a function of distance.  Its support remains within the support of R
// array<double, 2> energy(double dist, double R)
// // ----------------------------------------------------------------------------
// {
//   const double tmp = max(R-dist, double(0));
//   const double R3 = R*R*R; // normalizing factor so that E(0) = 1
//   return {tmp * tmp * tmp/R3, -3 * tmp*tmp/R3}; // energy and derivative
// }

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
