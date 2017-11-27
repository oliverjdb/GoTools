#include <iostream> /// @@ debug
#include <fstream> // @@ debug purposes
#include <stdexcept> // @@ debug
#include <set> // @@ debug
#include <algorithm>
#include <assert.h>
#include <stdexcept>
#include "triangulate_domain.h"
#include "tesselate_utils.h"
#include "basic_intersections.h"
#include "TriangleOctTree.h"

using namespace std;
using namespace TesselateUtils;

namespace {
  enum TriangleStatus {BEHIND_FRONT=0, PROBLEMATIC=1, FRONT_NON_DEL=2, FRONT_DEL=3};

  
  Triangle add_triangle(vector<Segment>& nsegs,   // input-output
			vector<Segment>& dsegs,   // input-output
			vector<uint>& unused_pts,    // input-output
			const Point2D* const points, // input only
			const double vdist);         // input only

  bool tet_found(const uint tri_ix,              
                 TriangleOctTree& tri_octtree,   // input-output
                 vector<TriangleStatus>& tri_flags,        // input-output
                 const vector<Triangle>& btris,
                 vector<uint>& unused_pts,
                 const Point3D* const points,
                 const uint num_bpoints,
                 const double vdist,
                 Tet& new_tet);                  // output
  
  void find_candidate_points(const Segment& s,
			     const vector<Segment>& nsegs,
			     const vector<uint>& unused_pts,
			     const Point2D* const points,
			     const double vdist,
			     vector<uint>& cand_pts,         // input-output
			     vector<uint>& all_neigh_pts);    // input-output

  void find_candidate_points(const Triangle& tri,
                             const TriangleOctTree& tri_octtree,
                             const vector<TriangleStatus>& tri_flags,
                             const vector<Triangle>& btris,                             
                             const vector<uint>& unused_pts,
                             const Point3D* const points,
                             const uint num_bpoints,
                             const double vdist,
                             vector<uint>& cand_pts,       // input-output
                             vector<uint>& all_neigh_pts,  // input-output
                             const bool is_non_del);
  uint best_candidate_point(const vector<uint>& cand_pts, // nb: will be modified
			    const Segment& seg,
			    const Point2D* const points);

  uint best_candidate_point(const vector<uint>& cand_pts, // nb: will be modified
                            const Triangle& cur_tri,
                            const Point3D* const points,
                            const bool is_obtuse_non_del);

  bool is_delaunay(const Triangle& tri,
		   const vector<uint>& neigh_pts,
		   const Point2D* const points);

  bool is_delaunay(const Triangle& tri, const uint chosen_pt,
                   const vector<uint>& neigh_pts, const Point3D* const points);
  
  bool add_or_remove_segment(const Segment& seg,
			     vector<Segment>& nsegs,
			     vector<Segment>& dsegs,
			     vector<uint>& unused_pts,
			     const bool delaunay,
			     const bool orientation); // 1: the new point is first in the segment
                                                      // 0: the new point is last in the segment

  bool add_or_remove_face(const Triangle& tri, // third corner of 'tri' is the new point
                          vector<Triangle>& ntris,
                          vector<Triangle>& dtris,
                          vector<Triangle>& ptris,
                          vector<uint>& unused_pts,
                          const bool delaunay);

  void add_face_to_front(const Triangle& tri,
                         TriangleOctTree& tri_octtree,
                         vector<TriangleStatus>& tri_flags,
                         const bool is_del);
  
  bool introduces_intersection(const Segment& s,
			       const uint pt_ix,
			       const vector<Segment>& nsegs,
			       const Point2D* const points,
			       const double tol);

  bool introduces_intersection(const Triangle& tri,
                               const uint pt_ix,
                               const TriangleOctTree& tri_octtree,
                               const vector<TriangleStatus>& tri_flags,
                               const bool is_non_del,
                               const Point3D* const points,
                               const double tol);

  bool search_and_erase_segment(const Segment& seg, vector<Segment>& segvec);
  bool search_and_erase_face(const Triangle& tri, vector<Triangle>& trivec);
  uint search_face(const Triangle& tri, const vector<Triangle>& trivec);

  uint max_index(const Triangle* const tris, uint num_tris);
  inline pair<Triangle, Triangle> boundary_triangles_of_edge(const Segment e,
                                                             const vector<Triangle>& btris,
                                                             uint& count);
  inline bool point_is_outside(const Point3D* const points, uint pt_ix, const Triangle& tri);
  array<Point3D, 3> modified_triangle(const array<Point3D, 3>& pts);
  uint purge_unused_triangles(TriangleOctTree& tri_octtree,
                              vector<TriangleStatus>& tri_flags);

};

namespace TesselateUtils {

// ============================================================================
// The algorithm used here is inspired by:
// S.H. Lo, "Delaunay Triangulation of Non-Convex Planar Domains" (1989)
vector<Triangle> triangulate_domain(const Point2D* const points,
				    const uint num_bpoints, 
				    const uint tot_num_points,
				    const double vdist)
// ============================================================================  
{

  // Preparing vectors keeping track of "non-delaunay" and delaunay segments.
  // The working front constists of both types of segments, but while adding
  // triangles we only need to check for intersections against the non-delaunay
  // segments.  At the beginning, the working front is equal to the boundary
  // polgygon, which is considered to consist solely of non-delaunay segments.
  // (If we were positive the boundary polygon was convex, we could probably
  // consider them delaunay segments instead).
  vector<Segment> nsegs(num_bpoints), dsegs; 

  // At start, all boundary segments are "non-delaunay".  We store them all.
  for (uint i = 0; i != num_bpoints; ++i)
    nsegs[i] = {i, (i+1)%num_bpoints};

  // We use this vector to keep track of indices to nodes that are still
  // interior or on the working front, which in the beginning means "all nodes".
  vector<uint> unused_pts(tot_num_points, 1);  // 1 for unused, 0 for used

  // Establish the result vector, and gradually fill it with rectangles until
  // there are no remaining points on the working front or in the interior.
  // NB: in the function call in the loop below, the three first arguments are
  // input-output (i.e. they will be modified by the function).
  vector<Triangle> result;

  while (nsegs.size() + dsegs.size() > 0)
    result.push_back(add_triangle(nsegs, dsegs, unused_pts, points, vdist));

  // sanity check: there should be no unused nodes left by now
  assert(accumulate(unused_pts.begin(), unused_pts.end(), 0) == 0);
  
  return result;
}

// ============================================================================  
void front_sanity_check(const vector<Triangle>& ntris,
                        const vector<Triangle>& dtris)
// ============================================================================  
{
  multiset<Segment> segs;
  auto extract_segments = [&] (const Triangle& t) {
    for (uint i = 0; i != 3; ++i) {
      uint a = min(t[i], t[(i+1)%3]);
      uint b = max(t[i], t[(i+1)%3]);
      segs.insert(Segment{a, b});
    }
  };
  for_each(ntris.begin(), ntris.end(), extract_segments);
  for_each(dtris.begin(), dtris.end(), extract_segments);

  for_each(segs.begin(), segs.end(), [&] (const Segment& s) {
      //assert(segs.count(s) == 2);
      if (segs.count(s) % 2 != 0) {
        //cout << "Something wronge happened.  Segs.count is: " << segs.count(s) << endl;
        //throw runtime_error("Something wrong happened.");
      }
    });
}
   
// ============================================================================
vector<Tet> construct_tets(const Point3D* const points,
                           const uint tot_num_points,
                           const Triangle* btris,
                           const uint num_btris,
                           const double vdist)
// ============================================================================
{
  TriangleOctTree tri_octtree(points, btris, num_btris);

  vector<TriangleStatus> tri_flags(num_btris, FRONT_NON_DEL); 
  
  // interior or on the working front, which in the beginning means "all nodes".
  vector<uint> unused_pts(tot_num_points, 1); // 1 for unused, 0 for used

  // determine how many boundary points there are 
  const uint num_bpoints = max_index(btris, num_btris) + 1;

  const vector<Triangle> btris_vec(btris, btris + num_btris);
  
  // Establish the result vector, and gradually fill it with tets until there
  // are no remaining points on the working front nor in the interior.
  vector<Tet> result;
  Tet new_tet;

  while (true) {
    // continue as long as there are unprocessed front triangles left

    // we seek to deplete non-delaunay triangles first, so choose one if there is one
    auto next = find(tri_flags.begin(), tri_flags.end(), FRONT_NON_DEL);
    if (next == tri_flags.end())
      // no non-delaunay triangles found.  Are there any delaunay triangles?
      next = find(tri_flags.begin(), tri_flags.end(), FRONT_DEL);
    
    if (next == tri_flags.end())
      // no more triangles left on the active front
      break;
    
    const uint ix = (uint) (next - tri_flags.begin());

    if (tet_found(ix, tri_octtree, tri_flags, btris_vec, unused_pts, points,
                  num_bpoints, vdist, new_tet))
      result.push_back(new_tet);

    cout << "Current number of tets: " << result.size() << '\n';
    cout << "Total number of triangles: " << tri_flags.size() << '\n';
    cout << "Non-del front triangles: " << count(tri_flags.begin(), tri_flags.end(), FRONT_NON_DEL) << '\n';
    cout << "Delaunay front triangles: " << count(tri_flags.begin(), tri_flags.end(), FRONT_DEL) << '\n';
    cout << "Bad triangles: " << count(tri_flags.begin(), tri_flags.end(), 1) << '\n';
    
  }

  if (count(tri_flags.begin(), tri_flags.end(), 2) > 0) {
    cout << "The following ptris were detected: \n";
    for (uint i = 0; i != (uint)tri_flags.size(); ++i) {
      if (tri_flags[i] == 2) {
        const Triangle& t = tri_octtree.getTri(i);
        cout << t[0] << " " << t[1] << " " << t[2] << "\n";
      }
    }
    cout << endl;
  }
  // sanity check: there should be no unused nodes left by now
  //  assert(accumulate(unused_pts.begin(), unused_pts.end(), 0) == 0);

  return result;
}

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
uint max_index(const Triangle* const tris, uint num_tris)
// ----------------------------------------------------------------------------
{
  uint result = 0;
  for (uint i = 0; i != num_tris; ++i) 
    result = max(result, *max_element(tris[i].begin(), tris[i].end()));
  return result;
}

// ----------------------------------------------------------------------------
bool tet_found(const uint tri_ix,
               TriangleOctTree& tri_octtree,
               vector<TriangleStatus>& tri_flags,
               const vector<Triangle>& btris,
               vector<uint>& unused_pts,
               const Point3D* const points,
               const uint num_bpoints,
               const double vdist,
               Tet& new_tet)
// ----------------------------------------------------------------------------
{
  // remove the triangle from the active front
  const TriangleStatus type = tri_flags[tri_ix];
  assert(type > 1);  // FRONT_NON_DEL (2) or FRONT_DEL (3)
  tri_flags[tri_ix] = BEHIND_FRONT;
  
  const Triangle cur_tri = tri_octtree.getTri(tri_ix);

  // If the triangle is non-delaunay and too obtuse, it may lead to the
  // generation of very bad tets.  In that case, we here choose triangle with a
  // similar non-obtuse one, and flag the resulting tet as non-delaunay.  Since
  // we are creating a non-delaunay tet, the intersection check must be done
  // against _all_ sufficiently close front triangles, not only non-delaunay
  // ones.  (We should not do this for slightly obtuse triangles - what we want
  // is to eliminate the worst cases where a circle centre is located far away
  // from the triangle corners in the triangle's plane, which may lead to
  // near-degenerate tets).
  bool is_obtuse_non_del = false;
  // const array<Point3D, 3> tripts {points[cur_tri[0]], points[cur_tri[1]], points[cur_tri[2]]};
  
  // if ((type == FRONT_NON_DEL) && is_obtuse_triangle(tripts[0], tripts[1], tripts[2])) {
                                                    
  //   // check if the triangle is too obtuse.  We here use the criteron that the
  //   // radius of the circumscribed cicle should not be larger than the longest
  //   // edge of the triangle.
  //   double r2;
  //   Point3D center;
  //   circumscribe_triangle(points[cur_tri[0]], points[cur_tri[1]], points[cur_tri[2]],
  //                         center, r2);
  //   is_obtuse_non_del = (r2 > 2 * norm2(tripts[1] - tripts[0]) ||
  //                        r2 > 2 * norm2(tripts[2] - tripts[1]) ||
  //                        r2 > 2 * norm2(tripts[0] - tripts[2]));
  // }
  
  // Identify all points within 'vdist' of triangle, and separate them into
  // 'candidate' and 'non-candidate' points, according to criteria analogous to
  // those mentioned in the comments to 'add_triangle'.  
  vector<uint> cand_pts, all_neigh_pts;
  find_candidate_points(cur_tri, tri_octtree, tri_flags, btris, unused_pts,
                        points, num_bpoints, vdist, cand_pts, all_neigh_pts,
                        type==FRONT_NON_DEL); 
  
  // The 'best' of the candidate points will be chosen as the fourth corner of
  // the tet.  The best point is the one with no other candidate points within
  // the circumscribing sphere of the generated tet.
  const uint chosen_pt = best_candidate_point(cand_pts, cur_tri, points, is_obtuse_non_del);

  // check if a candidate point was really found
  if (chosen_pt == (uint)cand_pts.size()) {
    // no point found.  No way to make a tet
    tri_flags[tri_ix] = PROBLEMATIC; // flag this triangle as problematic
    return false;
  }

  const bool is_del = (!is_obtuse_non_del) &&
                      is_delaunay(cur_tri, chosen_pt, all_neigh_pts, points);
    
  // We have defined the new tet.  Return it, and do all necessary
  // housekeeping
  new_tet = Tet {cur_tri[1], cur_tri[0], cur_tri[2], chosen_pt};
  
  // modify working front and remaining nodes.  Make sure to add them in
  // _clockwise_ corner order, since we want them to be pointing outwards of the
  // active front.
  add_face_to_front({cur_tri[0], cur_tri[1], chosen_pt}, tri_octtree, tri_flags, is_del);
  add_face_to_front({cur_tri[1], cur_tri[2], chosen_pt}, tri_octtree, tri_flags, is_del);
  add_face_to_front({cur_tri[2], cur_tri[0], chosen_pt}, tri_octtree, tri_flags, is_del);
  
  // Update which points remain on the active front.  Set all the involved
  // points to inactive, and add them back in if they are identified on the front again.
  unused_pts[chosen_pt] = 0;
  for (uint i = 0; i != 3; ++i) unused_pts[cur_tri[i]] = 0;
  for (uint i = 0; i != (uint)tri_flags.size(); ++i)
    if (tri_flags[i] != BEHIND_FRONT) {
      const Triangle& t = tri_octtree.getTri(i);
      for (uint c = 0; c != 3; ++c)
        unused_pts[t[c]] = 1;
    }

  // at periodic intervals, purge triangles that are no longer relevant
  // (i.e. those behind the front) from our data structures, to reduce search times
  const uint PURGE_INTERVAL = 500;
  if (tri_flags.size() % PURGE_INTERVAL == 0) {
    uint n = purge_unused_triangles(tri_octtree, tri_flags);
    cout << "Purged " << n  << " unused triangles." << endl;
    assert((uint)tri_flags.size() == tri_octtree.numTris());
  }
  
  return (chosen_pt < (uint)cand_pts.size());
}

// ----------------------------------------------------------------------------
uint purge_unused_triangles(TriangleOctTree& tri_octtree,
                            vector<TriangleStatus>& tri_flags)
// ----------------------------------------------------------------------------
{
  // identify indices of triangles no longer on the front
  vector<uint> ixs;
  vector<TriangleStatus> tri_flags_new;
  for (uint i = 0; i != (uint)tri_flags.size(); ++i)
    if (tri_flags[i] == BEHIND_FRONT)
      ixs.push_back(i);
    else
      tri_flags_new.push_back(tri_flags[i]);

  tri_flags.swap(tri_flags_new);

  // remove triangles from octtree
  tri_octtree.removeTris(ixs);
  
  // return the number of triangles removed
  return (uint)ixs.size();
}

  
// ----------------------------------------------------------------------------
uint search_face(const Triangle& tri, const vector<Triangle>& trivec)
// ----------------------------------------------------------------------------
{
  auto iter = find_if(trivec.begin(), trivec.end(), [&tri](Triangle t) {
      reverse(t.begin(), t.end());
      // flip triangle, rotate it three times, and check for a match each time
      for (uint i = 0; i!=3; ++i) {
        if (equal(t.begin(), t.end(), tri.begin()))
          return true;
        rotate(t.begin(), t.begin()+1,t.end());
      }
      
      return false;
    });
  return (uint) (iter - trivec.begin());
}
// ----------------------------------------------------------------------------
void add_face_to_front(const Triangle& tri,
                       TriangleOctTree& tri_octtree,
                       vector<TriangleStatus>& tri_flags,
                       const bool is_del)
// ----------------------------------------------------------------------------
{
  const uint tri_found = search_face(tri, *tri_octtree.getAllTris());
  if (tri_found < tri_flags.size()) {
    // triangle already exists.  Remove it from the front
    tri_flags[tri_found] = BEHIND_FRONT;
    return;
  }
  
  // This triangle has not yet been included.  Add it.
  tri_octtree.addTriangle(tri);
  tri_flags.push_back(is_del ? FRONT_DEL : FRONT_NON_DEL); // flag it as part of the active front
  
}

  
// ----------------------------------------------------------------------------
uint best_candidate_point(const vector<uint>& cand_pts,
                          const Triangle& cur_tri,
                          const Point3D* const points,
                          const bool is_obtuse_non_del)
// ----------------------------------------------------------------------------
{
  // find indices of candidate points
  vector<uint> cand_ixs;
  for (uint i = 0; i != (uint)cand_pts.size(); ++i)
    if (cand_pts[i])
      cand_ixs.push_back(i);

  if (cand_ixs.empty())
    // return an impossible index to signal that no candidate was found
    return (uint)cand_pts.size();

  // If the triangle is an obtuse, non-delaunay triangle, the delaunay tet may
  // be of very poor quality or degenerate.  In those case, replace the triangle
  // with a similar, right-angled one in the search below.  This will ensure
  // that centre of the circumscribed circle will not lie outside the triangle.
  array<Point3D,3> tripoints = {points[cur_tri[0]], points[cur_tri[1]], points[cur_tri[2]]};
  if (is_obtuse_non_del) {
    cout << "Triangle " << cur_tri << " is obtuse.   Replacing it." << endl;
    tripoints = modified_triangle(tripoints);
  }
    
  // searching for the candidate whose resulting tet has a minimal containing
  // sphere that does not contain any other candidate.
  Point3D center; // center of the minimal containing sphere
  double radius2; // squared radius of the containing sphere

  // @@ It may be faster to eliminate neighbors as one goes along (any neighbour
  // outside the sphere could be immediately disregarded and removed from the
  // vector).  This would however introduce a bit of extra logic.
  for (const auto ix : cand_ixs) {
    if (fitting_sphere(points[ix], tripoints[0], tripoints[1],
                       tripoints[2], center, radius2)) {
      // check whether the sphere contains any candidate point
      bool interior_point_found = false;
      for (const auto ix2 : cand_ixs)
        if (dist2(center, points[ix2]) < radius2) {
          interior_point_found = true;
          break;
        }
      if (!interior_point_found)
        return ix;
    }
  }
  // we should never get her by design, so throw an error if it happens
  throw runtime_error("No suitable best candidate.  This should not happen.");
  return -1; // dummy, to keep compiler happy
}

// ----------------------------------------------------------------------------  
void find_candidate_points(const Triangle& tri,
                           const TriangleOctTree& tri_octtree,
                           const vector<TriangleStatus>& tri_flags,
                           const vector<Triangle>& btris,
                           const vector<uint>& unused_pts,
                           const Point3D* const points,
                           const uint num_bpoints,
                           const double vdist,
                           vector<uint>& cand_pts,         // input-output
                           vector<uint>& all_neigh_pts,    // input-output
                           const bool is_non_del)
// ----------------------------------------------------------------------------
{
  const double TOL = 1.0e-10 * vdist; //1.0e-6; // @@ safe/general enough?
  const uint N = (uint)unused_pts.size();
  const Point3D& p1 = points[tri[0]]; // p1, p2 and p3 references are for 
  const Point3D& p2 = points[tri[1]]; // convenience only
  const Point3D& p3 = points[tri[2]];
  
  cand_pts.resize(N); fill(cand_pts.begin(), cand_pts.end(), 0);
  all_neigh_pts.resize(N); fill(all_neigh_pts.begin(), all_neigh_pts.end(), 0);
  
  for (uint i = 0; i != N; ++i) {
    if ((!unused_pts[i]) || (i== tri[0]) || (i == tri[1]) || (i == tri[2])) {
      // this is either a non-active point or one of the triangle corners; we
      // need not consider those.
      continue;
    }
    
    if (!point_on_triangle(points[i], p1, p2, p3, vdist))
      // point not close enough to triangle to be considered a neighbor
      continue;

    all_neigh_pts[i] = 1; // all points that are close enough, regardless of other factors
    if (i < num_bpoints) {
      // verify that we are not expanding the outer boundary by using this point
      // to construct a new tet
      bool ok = true;
      for (uint j = 0; j != 3; ++j) {
        if ((tri[j] < num_bpoints) && (tri[(j+1)%3] < num_bpoints))  {
          // this is a boundary edge, and we must ensure that the new point is
          // either already forming a triangle with this edge, or that any new
          // triangle formed will be in the interior of the object.
          uint count = 0;
          auto neigh_tris = boundary_triangles_of_edge(Segment {tri[j], tri[(j+1)%3]},
                                                       btris, count);
          if (count == 0)
            //no problem.  This is an internal edge
            continue;
          
          if ((find(neigh_tris.first.begin(), neigh_tris.first.end(), i) != neigh_tris.first.end()) ||
              (find(neigh_tris.second.begin(), neigh_tris.second.end(), i) != neigh_tris.second.end()))
            // No problem.  This triangle is already part of the boundary. 
            continue;
          
          // A new triangle will be introduced.  Verify that it will be an
          // internal triangle (i.e. not expand the existing boundary)
          if (point_is_outside(points, i, neigh_tris.first) &&
              point_is_outside(points, i, neigh_tris.second)) {
            ok = false;
            break;
          }
        }
      }
      if (!ok)
        // this point cannot be used.
        continue;
    }

    if ((point_on_inside_of_face(points[i], p1, p2, p3, TOL)) &&
        (!introduces_intersection(tri, i, tri_octtree, tri_flags,
                                  is_non_del, points, TOL)))
      cand_pts[i] = 1; // point is a candidate point
  }
}

// ----------------------------------------------------------------------------
inline bool point_is_outside(const Point3D* const points, uint pt_ix, const Triangle& tri)
// ----------------------------------------------------------------------------
{
  const Point3D v1 = points[tri[1]] - points[tri[0]];
  const Point3D v2 = points[tri[2]] - points[tri[0]];
  const Point3D v3 = points[pt_ix]  - points[tri[0]];
  return determinant3D(v1, v2, v3) >= 0;
}

// ----------------------------------------------------------------------------
inline pair<Triangle, Triangle> boundary_triangles_of_edge(const Segment e,
                                                           const vector<Triangle>& btris,
                                                           uint& count)
// ----------------------------------------------------------------------------
{
  count = 0;
  Triangle res[2];
  for (uint i = 0; i != (uint)btris.size(); ++i) 
    if ((e[0] == btris[i][0] || e[0] == btris[i][1] || e[0] == btris[i][2]) &&
        (e[1] == btris[i][0] || e[1] == btris[i][1] || e[1] == btris[i][2]))
      res[count++] = btris[i];

  assert((count == 0) || (count == 2));
  return pair<Triangle, Triangle> {res[0], res[1]};
}
  
// ----------------------------------------------------------------------------
bool introduces_intersection(const Triangle& tri,
                             const uint pt_ix,
                             const TriangleOctTree& tri_octtree,
                             const vector<TriangleStatus>& tri_flags,
                             const bool is_non_del,
                             const Point3D* const points,
                             const double tol)
// ----------------------------------------------------------------------------
{
  // Initial check: if there are no non-delaunay triangles on front, and if the
  // base triangle is also delaunay, there should be no possible intersections,
  // and we return directly
  if ((!is_non_del) &&
      find(tri_flags.begin(), tri_flags.end(), FRONT_NON_DEL) == tri_flags.end())
    return false;
  
  // checking each triangle for intersection against the triangles that would be
  // introduced if the point indiced by 'pt_ix' were to form a tet with the
  // already-existing triangle 'tri'.
  const array<array<Point3D, 3>, 3> new_tris {
    array<Point3D, 3> {points[tri[0]], points[tri[1]], points[pt_ix]},
    array<Point3D, 3> {points[tri[1]], points[tri[2]], points[pt_ix]},
    array<Point3D, 3> {points[tri[2]], points[tri[0]], points[pt_ix]}}; // introduced triangles

  vector<uint> candidates;

  for (uint i = 0; i != 3; ++i) { // loop over newly introduced triangles

    // check if the current new triangle will intersect any existing triangle on the front
    const array<Point3D, 3>& t = new_tris[i];
    tri_octtree.getIntersectionCandidates(t, candidates);

    vector<uint>::const_iterator cand_it;
    vector<TriangleStatus>::const_iterator flags_it;
    uint count = 0;
    for (flags_it = tri_flags.begin(), cand_it = candidates.begin(), count = 0;
         cand_it != candidates.end();
         ++flags_it, ++cand_it, ++count) {
      if ((*cand_it) && (((*flags_it) == FRONT_NON_DEL) ||
                         ((*flags_it) == PROBLEMATIC) ||
                         (is_non_del && (*flags_it == FRONT_DEL)))) {
        const Triangle& ft = tri_octtree.getTri(count); // triangle on existing front
        const array<Point3D, 3> tripoints {points[ft[0]], points[ft[1]], points[ft[2]]};
        if (isect_triangle_triangle_3D(&t[0], &tripoints[0], fabs(tol)) == OVERLAPPING)
          return true;
      }
    }
  }
  // no intersection detected
  return false;
}
    
// ----------------------------------------------------------------------------
Triangle add_triangle(vector<Segment>& nsegs,   
		      vector<Segment>& dsegs,   
		      vector<uint>& unused_pts,    
		      const Point2D* const points, 
		      const double vdist)
// ----------------------------------------------------------------------------
{
  // We choose the next segment to work with.  We aim to deplete the
  // non-delaunay segments as fast as possible, since we do not have to do
  // intersection checks againt the delaunay segments.  We therefore pick a
  // segment from 'nsegs' as long as there are any left.
  const Segment cur_seg = (nsegs.size() > 0) ? nsegs.back() : dsegs.back();
  (nsegs.size() > 0) ? nsegs.pop_back() : dsegs.pop_back();
  
  // Identify all points within 'vdist' of the segment.  We separate these into
  // 'candidate' and 'non-candidate' points.  Candidate points must fulfill the 
  // additional criteria:
  // (1) They have to be on the 'left side' of the current segment
  // (2) No intersections with other segments on on the working front must occur
  //     when connecting the candidate point to either endpoint of the current
  //     segment (sufficient to check agaisnt the non-delaunay segments)
  // We will chose the third point of the triangle from the candidate points.
  // We must still check against all neighbor points to determine whether
  // our triangle is delaunay or not (semi-delaunay is still "not delaunay")
  vector<uint> cand_pts, all_neigh_pts;
  find_candidate_points(cur_seg, nsegs, unused_pts, points, vdist,
			cand_pts, all_neigh_pts);

  // Choosing the best candidate point, i.e. a candidate points for which no
  // other candidate points are within the circumscribing circle of the triangle
  // generated if using this point.  
  const uint chosen_pt = best_candidate_point(cand_pts, cur_seg, points);

  // Determine whether the generated triangle is a "pure" delaunay triangle, or
  // if it contains other points caused by the nonconvexity of the working front
  // (or if it is semi-delaunay by having other points exactly on its
  // circumscription).
  const bool is_del = is_delaunay({cur_seg[0], cur_seg[1], chosen_pt},
				  all_neigh_pts, points);
  
  // Modify working front and remaining nodes
  const bool removed1 = add_or_remove_segment({cur_seg[0], chosen_pt}, nsegs, dsegs,
					      unused_pts, is_del, false);
  const bool removed2 = add_or_remove_segment({chosen_pt, cur_seg[1]}, nsegs, dsegs,
					      unused_pts, is_del, true);

  if (removed1 || removed2) {
    // at least one segment was removed. We update the list of active points
    unused_pts[cur_seg[0]] =  unused_pts[cur_seg[1]] = unused_pts[chosen_pt] = false;
    for (auto e : nsegs) 
      unused_pts[e[0]] = unused_pts[e[1]] = true;
    for (auto e : dsegs) 
      unused_pts[e[0]] = unused_pts[e[1]] = true;
  }
  
  return Triangle {cur_seg[0], cur_seg[1], chosen_pt};
}

// ----------------------------------------------------------------------------
void find_candidate_points(const Segment& s,
                           const vector<Segment>& nsegs,
                           const vector<uint>& unused_pts,
                           const Point2D* const points,
                           const double vdist,
                           vector<uint>& cand_pts,         // input-output
                           vector<uint>& all_neigh_pts)     // input-output
// ----------------------------------------------------------------------------
{
  const double TOL = 1.0e-9;  // @@ is this safe/general enough?
  const uint N = (uint)unused_pts.size(); // total number of points
  cand_pts.resize(N);     fill(cand_pts.begin(), cand_pts.end(), 0);
  all_neigh_pts.resize(N); fill(all_neigh_pts.begin(), all_neigh_pts.end(), 0);

  for(uint i = 0; i != N; ++i) {
    if ((!unused_pts[i]) || (i== s[0]) || (i == s[1])) {
      // this is either a non-active point or one of the segment endpoints; we
      // need not consider those
      continue;
    }
    if (point_on_line_segment(points[i], points[s[0]], points[s[1]], vdist, false)) {
      // point is close enough to be a neighbor point
      all_neigh_pts[i] = 1;
      // check if it is also a candidate point
      if ((projected_distance_to_line_2D(points[i], points[s[0]], points[s[1]]) < 0) &&
	  (!introduces_intersection(s, i, nsegs, points, TOL))) {
	cand_pts[i] = 1; // point is a candidate point
      } 
    }
  }
  
}

// ----------------------------------------------------------------------------  
bool introduces_intersection(const Segment& s,
			     const uint pt_ix,
			     const vector<Segment>& nsegs,
			     const Point2D* const points,
			     const double tol)
// ----------------------------------------------------------------------------  
{
  // @@ NB: there are algorithms for efficiently intersecting one set of line
  //    segments against another in one go.  If the below, "one-by-one"
  //    algorithm becomes a bottleneck, consider switching to the more efficient
  //    set-based approach.
  for (const auto& ns: nsegs)
    if ((ns[0] == pt_ix) || ns[1] == pt_ix)
      // candidate point is one of the endpoints of this segment.  No
      // new intersection possible.
      continue;
    else if ((ns[0] == s[0] && ns[1] == s[1]) || (ns[0] == s[1] && ns[1] == s[0]))
      // the two segments are the same - no new intersection created
      continue;
    else if (segments_intersect_2D(points[pt_ix], points[s[0]],
				   points[ns[0]], points[ns[1]], -tol) ||
	     segments_intersect_2D(points[pt_ix], points[s[1]],
				   points[ns[0]], points[ns[1]], -tol))
      return true;
  return false;
}


// ----------------------------------------------------------------------------  
uint best_candidate_point(const vector<uint>& cand_pts,
                          const Segment& seg,
                          const Point2D* const points)
// ----------------------------------------------------------------------------
{
  // find indices of candidate points
  vector<uint> cand_ixs;
  for (uint i = 0; i != (uint)cand_pts.size(); ++i)
    if (cand_pts[i])
      cand_ixs.push_back(i);

  if (cand_ixs.empty()) 
    // @@ this ought to be handled more gracefully
    throw runtime_error("No candidate found, increase 'vdist' and try again.");

  // Searching for the candidate whose resulting triangle has a circumscription
  // that does not contain any other candidate
  Point2D center;  // center of the circumscribed circle
  double radius2;   // squared radius of the circumscribed circle
  const Point2D& p1 = points[seg[0]];
  const Point2D& p2 = points[seg[1]];

  // avoid circles with radii large enough to cause numerical issues
  const double l2 = dist2(p1, p2); // squared length of segment
  
  // @@ It may be faster to eliminate neighbors as one goes along (any neighbour
  // outside the circumscribed circle could be disregarded and removed from the
  // vector).  This introduces a bit more logic, but something to consider if
  // this part of the code proves to be a bottleneck.
  for (const auto ix : cand_ixs)  {
    if (circumscribe_triangle(points[ix], p1, p2, center, radius2)) {
      // is the radius of the circle too big (points basically lying on a
      // straight line)?
      if (l2 < radius2 * numeric_limits<double>::epsilon())
        continue;  
      
      // check whether the circle contains any candidate point
      bool interior_point_found = false;
      for (const auto ix2 : cand_ixs) 
	if (dist2(center, points[ix2]) < radius2) {
	  interior_point_found = true;
	  break;
	}
      if (!interior_point_found)
	return ix;
    }
  }
  // we should never get here by design, so throw an error if it happens
  throw runtime_error("No suitable best candidate.  This should not happen.  "
		      "Please investigate.");
  return -1; // dummy, to keep compiler happy
}
  
// ----------------------------------------------------------------------------
bool is_delaunay(const Triangle& tri, 
                 const vector<uint>& neigh_pts,
                 const Point2D* const points)
// ----------------------------------------------------------------------------
{
  Point2D circ_center;
  double radius2;
  const double tol = 1e-2; //@@ OK?  In the worst case, we risk classifying a
			   //delaunay edge as non-delaunay, which would not do
			   //any harm to the final result, but slightly increase
			   //the computational cost (since there is now one more
			   //segment for which we need to check intersections).
  circumscribe_triangle(points[tri[0]], points[tri[1]], points[tri[2]],
			circ_center, radius2);

  radius2 *= (1+tol);  // slightly increase the radius to capture points on the
		      // edge of the circle.  If we capture points slightly
		      // outside, that should still be ok, see comment above.
  
  const uint N = (uint)neigh_pts.size();
  for (uint i = 0; i != N; ++i) 
    if ((neigh_pts[i]) &&                               // should be neigh point
	((i!=tri[0]) && (i!=tri[1]) && (i!=tri[2])) &&  // should not be corner of triangle
	(dist2(points[i], circ_center) < radius2))
      return false;
              
  return true;
}

// ----------------------------------------------------------------------------
bool is_delaunay(const Triangle& tri, const uint chosen_pt,
                 const vector<uint>& neigh_pts, const Point3D* const points)
// ----------------------------------------------------------------------------
{

  Point3D center;
  double radius2;
  const double tol = 1e-2; // same comment as for the 2D version of 'is_delaunay above
  fitting_sphere(points[tri[0]], points[tri[1]], points[tri[2]], points[chosen_pt],
                 center, radius2);
  radius2 *= (1+tol); // slightly increase radius to ensure points on the boundary are captured

  const uint N = (uint)neigh_pts.size();
  for (uint i = 0; i != N; ++i) 
    if ((neigh_pts[i]) && // should be a neighbor point
        (i != chosen_pt) && // should not be the chosen point
        ((i!=tri[0]) && (i!=tri[1]) && (i!=tri[2])) && // should not be part of triangle
        (dist2(points[i], center) <= radius2)) // should be within the radius
      return false;  // this point is inside the sphere, hence the proposed 'tet' is not delaunay

  return true;
}
  
// ----------------------------------------------------------------------------  
bool add_or_remove_segment(const Segment& seg,
                           vector<Segment>& nsegs,
                           vector<Segment>& dsegs,
                           vector<uint>& unused_pts,
                           const bool delaunay,
			   const bool orientation)
// ----------------------------------------------------------------------------
{
  // first, determine if segment already exists
  bool segment_found = search_and_erase_segment(seg, nsegs);
  if (!segment_found) {
    segment_found = search_and_erase_segment(seg, dsegs);
  }

  if (segment_found) {
    // segment already existed and has been removed (since both its triangles
    // have now been located).  We must remove the point that is no longer on
    // the working boundary.  @@ No!  That's not always the case.  There may be
    // other segments refering to this point if it connects two or more separate
    // 'holes' in the mesh.  Do not carry out the instruction below.

    //unused_pts[orientation ? seg[1] : seg[0]] = 0;
    return true;
  } else {
    // this is a newly introduced segment.  Add it to the respective list of
    // segments
    if (delaunay)
      dsegs.push_back(seg);
    else
      nsegs.push_back(seg);
  }
  return false;
}

// ----------------------------------------------------------------------------
bool search_and_erase_segment(const Segment& seg, vector<Segment>& segvec)
// ----------------------------------------------------------------------------
{
  auto iter = find_if(segvec.begin(), segvec.end(), [&seg](const Segment& s) {
      return (((s[0]==seg[0]) && (s[1]==seg[1])) ||
	      ((s[0]==seg[1]) && (s[1]==seg[0])));
    });

  if (iter != segvec.end()) {
    segvec.erase(iter);
    return true;
  }
  return false; 
}

// ----------------------------------------------------------------------------
bool add_or_remove_face(const Triangle& tri, // third corner of 'tri' is the new point
                        vector<Triangle>& ntris,
                        vector<Triangle>& dtris,
                        vector<Triangle>& ptris,
                        vector<uint>& unused_pts,
                        const bool delaunay)
// ----------------------------------------------------------------------------
{
  const bool tri_found = search_and_erase_face(tri, ntris) ||
                         search_and_erase_face(tri, dtris) ||
                         search_and_erase_face(tri, ptris);
  if (tri_found) 
    // triangle already existed and has been removed from the front.
    return true;

  // this is a newly introduced triangle.  Add it to the respective list of
  // triangles on the boundary
  (delaunay ? dtris : ntris).push_back(tri);
  return false;
}

// ----------------------------------------------------------------------------
bool search_and_erase_face(const Triangle& tri, vector<Triangle>& trivec)
// ----------------------------------------------------------------------------
{
  auto iter = find_if(trivec.begin(), trivec.end(), [&tri](Triangle t) {
      reverse(t.begin(), t.end());
      // flip triangle, rotate it three times, and check for a match each time
      for (uint i = 0; i!=3; ++i) {
        if (equal(t.begin(), t.end(), tri.begin()))
          return true;
        rotate(t.begin(), t.begin()+1,t.end());
      }
      
      return false;
    });
  if (iter != trivec.end()) {
    trivec.erase(iter);
    return true;
  }
  return false;
}

// ----------------------------------------------------------------------------
array<Point3D, 3> modified_triangle(const array<Point3D, 3>& pts)
// ----------------------------------------------------------------------------
{
  // compute a right-angled triangle to use instead of the obtuse triangle
  // defined by 'pts' when searching for a fourth point to construct a tet

  const array<Point3D, 3> v { pts[1] - pts[0], pts[2] - pts[1], pts[0] - pts[2]};
    
  const Point3D n = v[0] ^ v[1]; // non-normalized normal

  // compute longest edge
  const array<double, 3> l2 {norm2(v[0]), norm2(v[1]), norm2(v[2])};
  
  double cur_max = l2[0];
  int ix = 0;
  for (int i = 1; i != 3; ++i) {
    if (l2[i] > cur_max) {
      cur_max = l2[i];
      ix = i;
    }
  }
  // compute circle center and radius
  const double r = sqrt(cur_max)/2; // half the diameter
  const Point3D c = pts[ix] + 0.5 * v[ix];

  // defining first and second point of modified triangle
  array<Point3D, 3> result;
  result[0] = pts[ix];
  result[1] = pts[(ix+1)%3];
  
  // compute third point of modified triangle
  
  const Point3D radvec = result[1] - c;  // vector from circle center to second point
  result[2] = radvec ^ n; // cross product
  const double fac = r / norm(result[2]);
  result[2] *= fac;

  const Point3D rp = pts[(ix+2)%3]; // point to be removed
  if ( ( (rp-c) * result[2]) < 0) {
    result[2] *= -1;
    swap(result[0], result[1]);
  }
  
  result[2] += c;
  
  return result;
  
}
  
}; // end anonymous namespace

