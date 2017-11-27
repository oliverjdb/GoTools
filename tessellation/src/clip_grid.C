#include "clip_grid.h"
#include "tesselate_utils.h"

#include <fstream>  // @@ for debugging only
#include <iterator> // @@ for debugging only
using namespace std;
using namespace TesselateUtils;

namespace  {

void identify_intersected_cells(const array<double, 4>& bbox,
                                const array<uint, 2>& res,
                                const array<double, 2>& cell_len, 
                                const Point2D* const pcorners,
                                const uint num_pcorners,
                                ClippedDomainType* const result);

void identify_intersected_cells(const array<double, 6>& bbox,
                                const array<uint, 3>& res,
                                const array<double, 3>& cell_len, 
                                const Point3D* const pcorners,
                                const Triangle* const tris,
                                const uint num_tris,
                                ClippedDomainType* const result);

// 2D version
void classify_remaining_cells(const array<double, 4>& bbox,
                            const array<uint, 2>& res,
                            const Point2D* const pcorners,
                            const uint num_pcorners,
                            const double vdist,
                            ClippedDomainType* const result);

// 3D version
void classify_remaining_cells(const array<double, 6>& bbox,
                              const array<uint, 3>& res,
                              const Point3D* const pcorners,
                              const uint num_pcorners,
                              const Triangle* const tris,
                              const uint num_tris,
                              const double vdist,
                              ClippedDomainType* const result);
  
void identify_crossings(Point2D p1,
                        Point2D p2,
                        array<double, 4> bbox,
                        array<uint, 2> res,
                        array<double, 2> cell_len,
                        ClippedDomainType* const result,
                        const uint dir);


void reindex(vector<ClippedDomainType>& vec,
             const uint dim,
             const array<uint, 3> res);

array<Point2D, 2> tri_plane_isect_line(const array<Point3D, 3>& tricorners,
                                       const double plane_pval,
                                       const uint dim);
  
}; //end anonymous namespace

namespace TesselateUtils {

//----------------------------------------------------------------------------
ClippedGrid<2> clip_grid_polygon_2D(const Point2D* const pcorners,
                                    const uint num_pcorners,
                                    const double vdist, 
                                    const uint res_x,
                                    const uint res_y)
//----------------------------------------------------------------------------
{
  const auto bbox = bounding_box_2D(pcorners, num_pcorners);
  ClippedGrid<2> result { bbox,
                          {res_x, res_y},
                          {(bbox[1] - bbox[0])/res_x, (bbox[3] - bbox[2])/res_y},
                          vector<ClippedDomainType>(res_x * res_y, UNDETERMINED)};
  
  // identifying all cells intersected by polygon
  identify_intersected_cells(result.bbox, result.res, result.cell_len, pcorners,
                             num_pcorners, &(result.type[0]));

  // classify non-intersected cells as either 'OUTSIDE', 'FAR_INSIDE' or
  // 'CLOSE_INSIDE'
  classify_remaining_cells(result.bbox, result.res, pcorners,
                           num_pcorners, vdist, &(result.type[0]));
  
  return result;
};

// ----------------------------------------------------------------------------
ClippedGrid<3> clip_grid_shell_3D(const Point3D* const pcorners,
                                  const uint num_pcorners,
                                  const Triangle* const tris,
                                  const uint num_tris,
                                  const double vdist,
                                  const uint res_x,
                                  const uint res_y,
                                  const uint res_z)
// ----------------------------------------------------------------------------
{
  assert(res_x > 1 || res_y > 1 || res_z > 1); // necessary for this function to
                                               // work properly, since only
                                               // internal planes will be
                                               // checked for intersections
  
  const auto bbox = bounding_box_3D(pcorners, num_pcorners);
  //const auto bbox = bounding_box_3D(&krullcorners[0], num_pcorners);
  ClippedGrid<3> result {bbox,
                         {res_x, res_y, res_z},
                         {(bbox[1] - bbox[0])/ res_x,
                          (bbox[3] - bbox[2])/ res_y,
                          (bbox[5] - bbox[4])/ res_z},
                         vector<ClippedDomainType>(res_x * res_y * res_z, UNDETERMINED)};
    
  // identifying all cells intersected by a triangle
  identify_intersected_cells(result.bbox, result.res, result.cell_len, pcorners, 
                             tris, num_tris, &(result.type[0]));

  // Classify non-intersected cells as either 'OUTSIDE', 'FAR_INSIDE' or 'CLOSE_INSIDE'
  classify_remaining_cells(result.bbox, result.res, pcorners, num_pcorners, tris,
                           num_tris, vdist, &(result.type[0]));

  // ofstream os("krull"); //@@@
  // std::copy(result.type.begin(), result.type.end(), ostream_iterator<int>(os, " "));
  // os.close(); 
    
  
  return result;
}
  
}; //end namespace TesselateUtils

namespace {

//----------------------------------------------------------------------------  
void identify_intersected_cells(const array<double, 4>& bbox,
                                const array<uint, 2>& res,
                                const array<double, 2>& cell_len,
                                const Point2D* const pcorners,
                                const uint num_pcorners,
                                ClippedDomainType* const result)
//----------------------------------------------------------------------------  
{
  for (uint i = 0; i != num_pcorners; ++i) {
    const Point2D& p1 = pcorners[i];
    const Point2D& p2 = pcorners[(i+1) % num_pcorners];

    // identifying vertical and horizontal crossings
    identify_crossings(p1, p2, bbox, res, cell_len, result, 0); // vertical
    identify_crossings(p1, p2, bbox, res, cell_len, result, 1);  // horizontal
  }
}

//----------------------------------------------------------------------------
void identify_intersected_cells(const array<double, 6>& bbox,
                                const array<uint, 3>& res,
                                const array<double,3>& cell_len,
                                const Point3D* const pcorners,
                                const Triangle* const tris,
                                const uint num_tris,
                                ClippedDomainType* const result)
//----------------------------------------------------------------------------
{
  const double EPS = sqrt(numeric_limits<double>::epsilon());
  vector<ClippedDomainType> dresults(res[0] * res[1] * res[2], UNDETERMINED);
  for (uint dim = 0; dim != 3; ++dim) {
    // determine dimensionally reduced parameters
    const uint d1 = (dim+1)%3;
    const uint d2 = (dim+2)%3;
    const uint plane_size = res[d1] * res[d2]; 
    const double dlen = cell_len[dim];
    const uint dres = res[dim];
    const array<double, 4> dbox {bbox[2*d1], bbox[2*d1+1],
                                 bbox[2*d2], bbox[2*d2+1]};
    fill(dresults.begin(), dresults.end(), UNDETERMINED); // reset the temporary result vec.
    
    for (auto t = tris; t != tris + num_tris; ++t)  {

      // identify triangle corners
      const array<Point3D, 3> tc { pcorners[(*t)[0]],
                                   pcorners[(*t)[1]],
                                   pcorners[(*t)[2]] };

      // identify number of cells along direction 'dim' touched or intersected
      // by this triangle
      auto mm = minmax({tc[0][dim], tc[1][dim], tc[2][dim]}); 
      const double tol = EPS * (mm.second - mm.first);
      mm.first -= tol;
      mm.second += tol;
      const pair<uint, uint> range
        { min ((uint)floor(max(mm.first  - bbox[2*dim], 0.0)/dlen), dres - 1),
          min ((uint)floor(max(mm.second - bbox[2*dim], 0.0)/dlen), dres - 1) };

      for (uint pl_ix = range.first; pl_ix != range.second; ++pl_ix) {
        // we deal with the interface betweens cells in plane pl_ix and pl_ix+1
        const double plane_pval = bbox[2*dim] + (pl_ix+1) * dlen;
        const array<Point2D, 2> isect_seg = tri_plane_isect_line(tc, plane_pval, dim);

        // identifying the intersected cells in plane 'pl_ix'
        identify_intersected_cells(dbox, {res[d1], res[d2]}, {cell_len[d1], cell_len[d2]},
                                   &isect_seg[0], 2, &dresults[pl_ix * plane_size]);
        // we also need to include the same cells in plane 'pl_ix'+1
        identify_intersected_cells(dbox, {res[d1], res[d2]}, {cell_len[d1], cell_len[d2]},
                                   &isect_seg[0], 2, &dresults[(pl_ix+1) * plane_size]);
      }
    }
    
    // results are now stored as [d1, d2, dim], we need to re-shuffle results
    // according to original storage patterns before updating our final result vector
    reindex(dresults, dim, res);

    // transfer results to final result vector
    for (uint i = 0; i != (uint)dresults.size(); ++i)
      if (dresults[i] == INTERSECTED)
        result[i] = INTERSECTED;
  }
}

//----------------------------------------------------------------------------
array<Point2D, 2> tri_plane_isect_line(const array<Point3D, 3>& tc, // triangle corners
                                       const double plane_pval,
                                       const uint dim)
//----------------------------------------------------------------------------
{
  const double EPS = sqrt(numeric_limits<double>::epsilon());
  const auto bbox = bounding_box_3D(&tc[0], 3);
  const double TOL = EPS * max(1.0, fabs(plane_pval)); // @@ appropriate tolerance?
  
  // If two corners lie in the plane, return the segment between the two
  // corners.  Otherwise, if one corner lies in the plane, return a degenerated
  // segment consisting of that point twice.  (If this function was called, we
  // assume that not all three points are found in the plane).
  for (uint i = 0; i != 3; ++i) {
    if (fabs(tc[i][dim] - plane_pval) < TOL) {
      // is there another point touching the plane?
      for (uint j = i+1; j != 3; ++j) {
        if (fabs(tc[j][dim] - plane_pval) < TOL) {
          return array<Point2D, 2> { Point2D {tc[i][(dim+1)%3], tc[i][(dim+2)%3]},
                                     Point2D {tc[j][(dim+1)%3], tc[j][(dim+2)%3]} };
        }
      }
      return array<Point2D, 2> { Point2D {tc[i][(dim+1)%3], tc[i][(dim+2)%3]},
                                     Point2D {tc[i][(dim+1)%3], tc[i][(dim+2)%3]} };
    }
  }

  // If we got here, no corner lies in the plane, yet since this function was
  // called, there should be an intersection.  This means that there should be
  // exactly two sides of the triangle that intersect the plane.  We must
  // identify these, and return the points where they intersect the plane.
  uint num_found = 0;
  array<Point2D, 2> result;
  for (uint i = 0; i != 3; ++i) {
    const uint j = (i+1)%3;
    if ( (tc[i][dim] - plane_pval) * (tc[j][dim] - plane_pval) < 0 ) {
      
      const double s = (plane_pval - tc[i][dim]) / (tc[j][dim] - tc[i][dim]);
      const Point2D p1 {tc[i][(dim+1)%3], tc[i][(dim+2)%3]};
      const Point2D p2 {tc[j][(dim+1)%3], tc[j][(dim+2)%3]};
      result[num_found++] = p1 * (1-s) + p2 * s;
      if (num_found == 2)
        break;
    }
  }
  assert(num_found == 2);
  return result;
}
  
//----------------------------------------------------------------------------
void reindex(vector<ClippedDomainType>& vec, const uint dim, const array<uint, 3> res)
//----------------------------------------------------------------------------
{
  assert((uint)vec.size() == res[0] * res[1] * res[2]);
  const uint d1 = (dim + 1)%3;
  const uint d2 = (dim + 2)%3;
  vector<ClippedDomainType> tmp(vec.size());
  array<uint, 3> ijk;
  for (ijk[0] = 0; ijk[0] != res[0]; ++ijk[0])
    for (ijk[1] = 0; ijk[1] != res[1]; ++ijk[1])
      for (ijk[2] = 0; ijk[2] != res[2]; ++ijk[2])
        tmp[ijk[0] + res[0] * (ijk[1] + res[1] * ijk[2])] =
          vec[ijk[d1] + res[d1] * (ijk[d2] + res[d2] * ijk[dim])];

  swap(vec, tmp);
}
  
// //----------------------------------------------------------------------------
// void identify_cells_intersected_by_triangle(const array<double, 6>& bbox,
//                                             const array<uint, 3>& res,
//                                             const array<double, 3>& cell_len, 
//                                             const array<Point3D, 3>& t, // triangle
//                                             vector<ClippedDomainType>& result)
// //----------------------------------------------------------------------------
// {
//   const double EPS = sqrt(numeric_limits<double>::epsilon());
//   for (uint dim = 0; dim != 3; ++dim) {
//     const double len  = cell_len[dim];
//     // identify the number of planes along direction 'dim' that are touched or
//     // intersected by the triangle
//     auto mm = minmax({t[0][dim], t[1][dim], t[2][dim]}); 
//     const double tol = EPS * (mm.second - mm.first);
//     mm.first -= tol;
//     mm.second += tol;
//     const pair<uint, uint> range
//     { min ((uint)floor(max(mm.first  - bbox[dim], 0.0)/len), res[dim] - 1),
//       min ((uint)floor(max(mm.second - bbox[dim], 0.0)/len), res[dim] - 1) };

//     // loop across each identified plane
//     for (uint pl_ix = range.first; pl_ix != range.second; ++pl_ix) {

//       const double plane_pval = (pl_ix+1) * len;
//       array<Point2D, 2> isect_seg = tri_plane_isect_line(t, plane_pval, dim);
      
//       // solve lower-dimensional problem
      
//       identify_intersected_cells(
      
//     }
//   }
// }
  
//----------------------------------------------------------------------------
void classify_remaining_cells(const array<double, 4>& bbox,
                              const array<uint, 2>& res,
                              const Point2D* const pcorners,
                              const uint num_pcorners,
                              const double vdist,
                              ClippedDomainType* const result)
//----------------------------------------------------------------------------
{
  const double dx = (bbox[1] - bbox[0]) / (double)res[0];
  const double dy = (bbox[3] - bbox[2]) / (double)res[1];
  for (uint y_ix = 0; y_ix != res[1]; ++y_ix) {
    for (uint x_ix = 0; x_ix != res[0]; ++x_ix) {
      const uint ix = y_ix * res[0] + x_ix;
      if (result[ix] == UNDETERMINED) {
        // This cell is not intersected by the polygon.  Check if it falls outside.
        const Point2D centroid { bbox[0] + (0.5 + (double)x_ix) * dx,
                                 bbox[2] + (0.5 + (double)y_ix) * dy};
        uint dummy = 0; //we do not need this value, but required for function call below
        const double r = max(dx, dy);
        result[ix] =
          (!inpolygon(centroid, &pcorners[0], num_pcorners, 0.0, (bool&)dummy)) ? OUTSIDE:
          ((dist2(centroid,
                (closest_point_on_loop(centroid,
                                              pcorners,
                                              num_pcorners,
                                              dummy))) + r*r) < vdist*vdist) ? CLOSE_INSIDE :
                                                                               FAR_INSIDE;
      }
    }
  }
}

//----------------------------------------------------------------------------
void classify_remaining_cells(const array<double, 6>& bbox,
                              const array<uint, 3>& res,
                              const Point3D* const pcorners,
                              const uint num_pcorners,
                              const Triangle* const tris,
                              const uint num_tris,
                              const double vdist,
                              ClippedDomainType* const result)
//----------------------------------------------------------------------------
{
  const double dx = (bbox[1] - bbox[0]) / (double)res[0];
  const double dy = (bbox[3] - bbox[2]) / (double)res[1];
  const double dz = (bbox[5] - bbox[4]) / (double)res[2];  
  const uint num_cells = res[0] * res[1] * res[2];
  for (uint ix = 0; ix != num_cells; ++ix) {
    if (result[ix] == UNDETERMINED) {
      // This cell is not intersected by any triangle.  Check whether it falls
      // outside
      const uint z_ix = (uint)floor(ix / (res[0] * res[1]));
      const uint tmp = ix % (res[0] * res[1]);
      const uint y_ix = (uint)floor(tmp / res[0]);
      const uint x_ix = tmp % res[0];
      const Point3D centroid { bbox[0] + (0.5 + (double)x_ix) * dx,
                               bbox[2] + (0.5 + (double)y_ix) * dy,
                               bbox[4] + (0.5 + (double)z_ix) * dz};
      uint dummy = 0; // value not needed, but needed for function call below
      const double r = max({dx, dy, dz});
      result[ix] =
        (!inside_shell(centroid, pcorners, tris, num_tris, 0.0, (bool&)dummy)) ? OUTSIDE:
        ((dist2(centroid,
                (closest_point_on_triangle_surface(centroid,
                                                   pcorners,
                                                   num_pcorners,
                                                   tris,
                                                   num_tris,
                                                   dummy))) + r*r) < vdist*vdist) ? CLOSE_INSIDE:
                                                                                    FAR_INSIDE;
    }
  }
}
                               
//----------------------------------------------------------------------------  
void identify_crossings(Point2D p1,
                        Point2D p2,
                        array<double, 4> bbox,
                        array<uint, 2> res,
                        array<double, 2> cell_len,
                        ClippedDomainType* const result,
                        const uint dir)
//----------------------------------------------------------------------------
{
  assert(dir==0 || dir == 1);
  // ensure the direction we will examine is the first one
  if (dir==1) {
    swap(p1[0], p1[1]);
    swap(p2[0], p2[1]);
    swap(bbox[0], bbox[2]); swap(bbox[1], bbox[3]);
    swap(res[0], res[1]);
    swap(cell_len[0], cell_len[1]);
  }

  // const array<double, 2> cell_len = {(bbox[1] - bbox[0])/res[0],
  //                                    (bbox[3] - bbox[2])/res[1]};

  // Determine the range of cells along the chosen direction that are
  // intersected by the line segment from p1 to p2
  auto mm = minmax({p1[0], p2[0]});
  const double EPS = sqrt(numeric_limits<double>::epsilon());
  // add small tolerance to avoid missing edge cases
  const double tol = EPS * (mm.second - mm.first);
  mm.first -= tol;
  mm.second += tol;
  const pair<uint, uint> ix_range
     {min((uint)floor(max(mm.first  - bbox[0], 0.0)/cell_len[0]), res[0] - 1),
      min((uint)floor(max(mm.second - bbox[0], 0.0)/cell_len[0]), res[0] - 1)};
  
  vector<array<uint,2>> ixs_pairs;
  for (uint i = ix_range.first; i != ix_range.second; ++i) {
    const double pval = (i+1) * cell_len[0] + bbox[0];

    // determine at what point, in the other direction than 'dir', does the
    // intersection occur.
    const double t = (pval - p2[0]) / (p1[0] - p2[0]);
    
    const double pval_other = t * p1[1] + (1-t) * p2[1];

    // determine the second index of the cell
    const double pos = (pval_other - bbox[2])/cell_len[1]; // @@ won't work for cell_len[1] = 0
    const uint other_ix = min((uint)max(floor(pos), 0.0), res[1] - 1);
    
    ixs_pairs.push_back({i, other_ix});

    // handling edge cases
    if (floor(pos-EPS) < other_ix && other_ix > 0) {
      ixs_pairs.push_back( {i, min((uint)max(floor(pos-EPS), 0.0), res[1] - 1)});
    } else if (floor(pos+EPS) > other_ix && other_ix < res[1] - 1) {
      ixs_pairs.push_back( {i, min((uint)max(floor(pos+EPS), 0.0), res[1] - 1)});
    }
  }

  for (uint i = 0; i != ixs_pairs.size(); ++i) {
    array<uint, 2> ix_pair = ixs_pairs[i]; 
    uint ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];
    
    result[ix] = INTERSECTED;
    
    ix_pair[0] += 1;
    ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];

    result[ix] = INTERSECTED;
  }
}

  
}; // end anonymous namespace 
