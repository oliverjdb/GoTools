#ifndef _INTERPOINT_DISTANCES_IMPL_H
#define _INTERPOINT_DISTANCES_IMPL_H

#include <tuple>
#include <stdexcept>
#include "tesselate_utils.h"
#include "interpoint_distances.h"

namespace {
template<typename PointXD> std::vector<TesselateUtils::DistanceEntry> 
interpoint_distances_bruteforce_impl(const PointXD* const points,
				     const uint num_points,
				     const double R);

template<typename PointXD> std::vector<TesselateUtils::DistanceEntry> 
interpoint_distances_smart_impl(const PointXD* const points,
				const uint num_points,
				const double R);

template<typename PointXD>  
struct BinnedPoints {
  std::vector<PointXD> binned_points;
  std::vector<uint> original_index;
  std::vector<uint> indices;
  uint num_bins_x;
  uint num_bins_y;
  uint num_bins_z; // only used for PointXD = Point3D
};

template<typename PointXD>
BinnedPoints<PointXD>  bin_points(const PointXD* const points,
                                  const uint num_points,
                                  const double R);

void add_distances_from_bins(const BinnedPoints<TesselateUtils::Point2D>& binned_points,
			     const uint ix1, const uint iy1,
			     const uint ix2, const uint iy2,
			     const double R, 
			     std::vector<TesselateUtils::DistanceEntry>& result);
void add_distances_from_bins(const BinnedPoints<TesselateUtils::Point3D>& binned_points,
			     const uint ix1, const uint iy1, const uint iz1,
			     const uint ix2, const uint iy2, const uint iz2,
                             const double R, 
			     std::vector<TesselateUtils::DistanceEntry>& result);

};


namespace TesselateUtils {


// ============================================================================
template<typename PointXD>  
std::vector<DistanceEntry> interpoint_distances(const PointXD* const points,
                                                const uint num_points,
                                                const double R,
                                                const bool bruteforce)
// ============================================================================
{
  if (bruteforce) {
    // O(N^2) complexity, but competitive or even faster on small data
    return interpoint_distances_bruteforce_impl(points, num_points, R);
  }
  // O(N log N) complexity
  return interpoint_distances_smart_impl(points, num_points, R);
}
  
// ============================================================================
template<typename PointXD>
std::vector<DistanceEntry> interpoint_distances(const PointXD* points,
                                                const uint num_points,
                                                const PointXD& p,
                                                const double R)
// ============================================================================
{
  const double R2 = R*R; 
  std::vector<DistanceEntry> result;
  const auto points_end = points + num_points;
  for (auto q = points; q != points_end; ++q)
    if (dist2(p, *q) < R2)
      result.push_back( {uint(q-points), 0, dist(p, *q)});
  return result;
}

// ============================================================================  
template<typename PointXD>
std::vector<DistanceEntry> interpoint_distances(const PointXD* points,
                                                const uint num_points,
                                                const PointXD* other_points,
                                                const uint num_other_points,
                                                const double R)
// ============================================================================  
{
  const double R2 = R*R;
  std::vector<DistanceEntry> result;
  result.reserve(1000);  // hopefully sufficient to avoid std::vector reallocation
  
  // We only use a bruteforce implementation here, since we do not expect
  // 'num_points' to be particularly big in this case. 
  const auto points_end = points + num_points;
  const auto other_points_end = other_points + num_other_points;
  
  for (auto p = points; p != points_end; ++p)
    for (auto q = other_points; q != other_points_end; ++q)
      if ((dist2(*p, *q) < R2) && p!=q) // skip if both points are the same
	result.push_back({uint(p-points), uint(q-other_points), dist(*p, *q)});

  return result;
}

}; // end namespace TesselateUtils




namespace {

// ----------------------------------------------------------------------------
template<> BinnedPoints<TesselateUtils::Point2D> 
bin_points(const TesselateUtils::Point2D* const points,
           const uint num_points,
           const double R)
// ----------------------------------------------------------------------------
{
  const auto bbox = TesselateUtils::bounding_box_2D(points, num_points); 
  const double xmin = bbox[0]; const double xmax = bbox[1];
  const double ymin = bbox[2]; const double ymax = bbox[3];

  std::vector<std::array<uint, 2>> bins(num_points); // identify the bin of each point
  const double Rinv = 1.0/R;
  for (uint i = 0; i != num_points; ++i) {
    bins[i][0] = (uint)ceil((points[i][0] - xmin) * Rinv);
    bins[i][1] = (uint)ceil((points[i][1] - ymin) * Rinv);
  }
  const uint num_bins_x = (uint)ceil((xmax-xmin) * Rinv) + 1;
  const uint num_bins_y = (uint)ceil((ymax-ymin) * Rinv) + 1;

  // Sort points according to bin
  BinnedPoints<TesselateUtils::Point2D>
    result { {num_points, {0, 0}}, std::vector<uint>(num_points, 0),
             {}, num_bins_x, num_bins_y, 1};
  uint pos = 0;
  result.indices.push_back(pos);
  
  for (uint iy = 0; iy != num_bins_y; ++iy) {
    for (uint ix = 0; ix != num_bins_x; ++ix) {
      for (auto& b : bins) {
	if ((b[0] == ix) && (b[1] == iy)) {
	  const uint orig_ix = (uint)(&b-&bins[0]);
	  result.original_index[pos] = orig_ix;
	  result.binned_points[pos++] = points[orig_ix];
	}
      }
      result.indices.push_back(pos);
    }
  }
  return result;
}
    
// ----------------------------------------------------------------------------
template<> BinnedPoints<TesselateUtils::Point3D>
bin_points(const TesselateUtils::Point3D* const points,
           const uint num_points,
           const double R)
// ----------------------------------------------------------------------------
{
  const auto bbox = TesselateUtils::bounding_box_3D(points, num_points); 
  const double xmin = bbox[0]; const double xmax = bbox[1];
  const double ymin = bbox[2]; const double ymax = bbox[3];
  const double zmin = bbox[4]; const double zmax = bbox[5];

  std::vector<std::array<uint, 3>> bins(num_points); // identify the bin of each point
  const double Rinv = 1.0/R;
  for (uint i = 0; i != num_points; ++i) {
    bins[i][0] = (uint)ceil((points[i][0] - xmin) * Rinv);
    bins[i][1] = (uint)ceil((points[i][1] - ymin) * Rinv);
    bins[i][2] = (uint)ceil((points[i][2] - zmin) * Rinv);
  }

  const uint num_bins_x = (uint)ceil((xmax-xmin) * Rinv) + 1;
  const uint num_bins_y = (uint)ceil((ymax-ymin) * Rinv) + 1;
  const uint num_bins_z = (uint)ceil((zmax-zmin) * Rinv) + 1;
  
  // Sort points according to bin
  BinnedPoints<TesselateUtils::Point3D> result {
    {num_points, {0, 0}},
    std::vector<uint>(num_points, 0),
    {},
    num_bins_x, num_bins_y, num_bins_z
  };
  
  uint pos = 0;
  result.indices.push_back(pos);

  for (uint iz = 0; iz != num_bins_z; ++iz) {
    for (uint iy = 0; iy != num_bins_y; ++iy) {
      for (uint ix = 0; ix != num_bins_x; ++ix) {
        for (auto& b : bins) {
          if ((b[0] == ix) && (b[1] == iy) && (b[2] == iz)) {
            const uint orig_ix = (uint)(&b-&bins[0]);
            result.original_index[pos] = orig_ix;
            result.binned_points[pos++] = points[orig_ix];
          }
        }
        result.indices.push_back(pos);
      }
    }
  }
  return result;
}
    
// ----------------------------------------------------------------------------  
template<typename PointXD> std::vector<TesselateUtils::DistanceEntry> 
interpoint_distances_bruteforce_impl(const PointXD* const points, 
				     const uint num_points,
				     const double R)
// ----------------------------------------------------------------------------
{
  const double R2 = R*R;
  std::vector<TesselateUtils::DistanceEntry> result;
  // brute force implementation
  // @@ This is an N^2 algorithm, so is not expected to scale wells

  const auto points_end = points + num_points;
  for (auto p = points; p != points_end; ++p)
    for (auto q = p+1; q != points_end; ++q)
      if (TesselateUtils::dist2(*p, *q) < R2)
  	result.push_back({uint(p-points),
                          uint(q-points),
                          TesselateUtils::dist(*p, *q)});
  return result;  
}

// ----------------------------------------------------------------------------
template<>  // specialization for Point2D
std::vector<TesselateUtils::DistanceEntry> 
interpoint_distances_smart_impl(const TesselateUtils::Point2D* const points, 
				const uint num_points, 
				const double R)
// ----------------------------------------------------------------------------
{
  // sort points in bins depending on spatial position
  const BinnedPoints<TesselateUtils::Point2D> binned_points = bin_points(points, num_points, R);
  
  std::vector<TesselateUtils::DistanceEntry> result;  
  const uint num_bins_x = binned_points.num_bins_x;
  const uint num_bins_y = binned_points.num_bins_y;
  
  // testing bin against themselves and neighbors
  for (uint iy = 0; iy != num_bins_y; ++iy) {
    for (uint ix = 0; ix != num_bins_x; ++ix) {
      // comparing current bin againts itself
      add_distances_from_bins(binned_points, ix, iy, ix, iy, R, result);
      if (ix < num_bins_x-1) // compare against right neighbor
	add_distances_from_bins(binned_points, ix, iy, ix+1, iy, R, result);
      if (iy < num_bins_y-1) // compare against top neighbor
	add_distances_from_bins(binned_points, ix, iy, ix, iy+1, R, result);
      if ((ix < num_bins_x-1) && (iy < num_bins_y-1)) // diagonal compare 1
	add_distances_from_bins(binned_points, ix, iy, ix+1, iy+1, R, result);
      if ((ix > 0) && (iy < num_bins_y-1)) // diagonal compare 2
	add_distances_from_bins(binned_points, ix, iy, ix-1, iy+1, R, result);
    }
  }
  return result;
}


// ----------------------------------------------------------------------------
template<>  // specialization for TesselateUtils::Point3D
std::vector<TesselateUtils::DistanceEntry> 
interpoint_distances_smart_impl(const TesselateUtils::Point3D* const points, 
				const uint num_points, 
				const double R)
// ----------------------------------------------------------------------------
{
  // sort points in bins depending on spatial position
  const BinnedPoints<TesselateUtils::Point3D> binned_points = bin_points(points, num_points, R);
  
  std::vector<TesselateUtils::DistanceEntry> result;  
  const uint num_bins_x = binned_points.num_bins_x;
  const uint num_bins_y = binned_points.num_bins_y;
  const uint num_bins_z = binned_points.num_bins_z;
  
  // testing bin against themselves and neighbors
  for (uint iz = 0; iz != num_bins_z; ++iz) {
    for (uint iy = 0; iy != num_bins_y; ++iy) {
      for (uint ix = 0; ix != num_bins_x; ++ix) {
        // comparing current bin againts itself
        add_distances_from_bins(binned_points, ix, iy, iz, ix, iy, iz, R, result);
        if (ix < num_bins_x-1) // compare against right neighbor
          add_distances_from_bins(binned_points, ix, iy, iz, ix+1, iy, iz, R, result);
        if (iy < num_bins_y-1) // compare against top neighbor
          add_distances_from_bins(binned_points, ix, iy, iz, ix, iy+1, iz, R, result);
        if ((ix < num_bins_x-1) && (iy < num_bins_y-1)) // diagonal compare 1
          add_distances_from_bins(binned_points, ix, iy, iz, ix+1, iy+1, iz, R, result);
        if ((ix > 0) && (iy < num_bins_y-1)) // diagonal compare 2
          add_distances_from_bins(binned_points, ix, iy, iz, ix-1, iy+1, iz, R, result);

        // compares against bins in the next z-layer
        if (iz != num_bins_z - 1) {
          const std::array<int, 2> xinc_range {(ix > 0) ? -1 : 0, (ix < num_bins_x-1) ? 1 : 0};
          const std::array<int, 2> yinc_range {(iy > 0) ? -1 : 0, (iy < num_bins_y-1) ? 1 : 0};
            
          for (int xinc = xinc_range[0]; xinc != xinc_range[1]+1; ++xinc)
            for (int yinc = yinc_range[0]; yinc != yinc_range[1]+1; ++yinc)
              add_distances_from_bins(binned_points, ix, iy, iz,
                                      ix + xinc, iy + yinc, iz + 1, R, result);
        }
      }
    }
  }
  return result;
}

// ----------------------------------------------------------------------------
inline void add_distances_from_bins(const BinnedPoints<TesselateUtils::Point2D>& binned_points,
                                    const uint ix1, const uint iy1,
                                    const uint ix2, const uint iy2,
                                    const double R, 
                                    std::vector<TesselateUtils::DistanceEntry>& result)
// ----------------------------------------------------------------------------  
{
  const uint bin_1_linear_ix = iy1 * binned_points.num_bins_x + ix1; 
  const uint bin_2_linear_ix = iy2 * binned_points.num_bins_x + ix2;

  // if we are comparing a bin against itself, we must take care not to register
  // each distance twice
  const bool same_bin = bin_1_linear_ix == bin_2_linear_ix;

  const uint start_ix_1 = binned_points.indices[bin_1_linear_ix];
  const uint num_pts_1 = binned_points.indices[bin_1_linear_ix + 1] - start_ix_1;
  const uint start_ix_2 = binned_points.indices[bin_2_linear_ix];;
  const uint num_pts_2 = binned_points.indices[bin_2_linear_ix + 1] - start_ix_2;

  auto new_dists = (same_bin) ?
    interpoint_distances(&binned_points.binned_points[start_ix_1], num_pts_1, R, true) :
    interpoint_distances(&binned_points.binned_points[start_ix_1], num_pts_1,
			 &binned_points.binned_points[start_ix_2], num_pts_2, R);

  // indices in new_dists refer to the indices relative to the local bin.  We
  // here convert to global indices
  for (auto& d : new_dists) {
    d.p1_ix = binned_points.original_index[start_ix_1 + d.p1_ix];
    d.p2_ix = binned_points.original_index[start_ix_2 + d.p2_ix];
  }

  result.insert(result.end(), new_dists.begin(), new_dists.end());  
}

// ----------------------------------------------------------------------------
inline void add_distances_from_bins(const BinnedPoints<TesselateUtils::Point3D>& binned_points,
                                    const uint ix1, const uint iy1, const uint iz1,
                                    const uint ix2, const uint iy2, const uint iz2,
                                    const double R, 
                                    std::vector<TesselateUtils::DistanceEntry>& result)
// ----------------------------------------------------------------------------  
{
  const uint bin_1_linear_ix = (iz1 * binned_points.num_bins_y + iy1) * binned_points.num_bins_x + ix1; 
  const uint bin_2_linear_ix = (iz2 * binned_points.num_bins_y + iy2) * binned_points.num_bins_x + ix2;

  // if we are comparing a bin against itself, we must take care not to register
  // each distance twice
  const bool same_bin = bin_1_linear_ix == bin_2_linear_ix;

  const uint start_ix_1 = binned_points.indices[bin_1_linear_ix];
  const uint num_pts_1 = binned_points.indices[bin_1_linear_ix + 1] - start_ix_1;
  const uint start_ix_2 = binned_points.indices[bin_2_linear_ix];;
  const uint num_pts_2 = binned_points.indices[bin_2_linear_ix + 1] - start_ix_2;

  auto new_dists = (same_bin) ?
    interpoint_distances(&binned_points.binned_points[start_ix_1], num_pts_1, R, true) :
    interpoint_distances(&binned_points.binned_points[start_ix_1], num_pts_1,
			 &binned_points.binned_points[start_ix_2], num_pts_2, R);

  // indices in new_dists refer to the indices relative to the local bin.  We
  // here convert to global indices
  for (auto& d : new_dists) {
    d.p1_ix = binned_points.original_index[start_ix_1 + d.p1_ix];
    d.p2_ix = binned_points.original_index[start_ix_2 + d.p2_ix];
  }

  result.insert(result.end(), new_dists.begin(), new_dists.end());  
}

};



#endif
