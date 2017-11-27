#include "TriangleOctTree.h"
#include <assert.h>

using namespace std;
using namespace TesselateUtils;

namespace {

// ----------------------------------------------------------------------------
inline uint max_ix(const Triangle* const tris, uint num_tris)
// ----------------------------------------------------------------------------
{

  const auto me = max_element(tris, tris + num_tris,
                        [](const Triangle& t1, const Triangle& t2) {
                           return max_element(t1.begin(), t1.end()) <
                                  max_element(t2.begin(), t2.end());});
  return *max_element(me->begin(), me->end());
}

// ----------------------------------------------------------------------------
inline array<double, 3> compute_midvals(const array<double, 6>& bbox)
// ----------------------------------------------------------------------------
{
  array<double, 3> result;
  for (uint i = 0; i != 3; ++i)
    result[i] = 0.5 * (bbox[2*i] + bbox[2*i+1]);
  return result;
}

// ----------------------------------------------------------------------------  
inline array<double, 6> sub_box(const array<double, 6>& bbox, uint oct_ix)
// ----------------------------------------------------------------------------  
{ 
  array<double, 6> result;

  // xmin or xmax?
  if (oct_ix%2 == 0) {
    result[0] = bbox[0];
    result[1] = 0.5 * (bbox[0] + bbox[1]);
  } else { 
    result[0] = 0.5 * (bbox[0] + bbox[1]);
    result[1] = bbox[1];
  }

  // ymin or ymax?
  if ( (oct_ix/2) % 2 == 0) {
    result[2] = bbox[2];
    result[3] = 0.5 * (bbox[2] + bbox[3]);
  } else {
    result[2] = 0.5 * (bbox[2] + bbox[3]);
    result[3] = bbox[3];
  }
  
  // zmin or zmax?
  if ( (oct_ix / 4) == 0) {
    result[4] = bbox[4];
    result[5] = 0.5 * (bbox[4] + bbox[5]);
  } else {
    result[4] = 0.5 * (bbox[4] + bbox[5]);
    result[5] = bbox[5];
  }
  return result;
}
  
}; // end anonymous namespace

namespace TesselateUtils
{

// ============================================================================
TriangleOctTree::TriangleOctTree(const Point3D* const pts,
                                 const Triangle* const tris,
                                 const uint num_tris,
                                 const uint max_num_ixs,
                                 const uint max_depth)
// ============================================================================
  : max_num_ixs_(max_num_ixs), // how many indices the volume can contain before it subdivides
    cur_depth_(0),
    max_depth_(max_depth),
    points_(pts),
    tris_(shared_ptr<vector<Triangle>>(new vector<Triangle>(tris, tris + num_tris))),
    bbox_(bounding_box_3D(pts, max_ix(tris, num_tris))),
    midvals_(compute_midvals(bbox_)), 
    indices_(num_tris)
{
  for (uint i = 0; i != num_tris; ++i)
    indices_[i] = i; // all triangles belong to this oct-tree

  reorganize_if_necessary();
}

// ============================================================================
void TriangleOctTree::addTriangle(const Triangle& t)
// ============================================================================
{
  tris_->push_back(t); // add the triangle itself
  include_last_triangle();  // update the indices_ (of itself or of its children)
  reorganize_if_necessary();
}

// ============================================================================
void TriangleOctTree::removeTris(const vector<uint>& ixs)
// ============================================================================
{
  if (is_subdivided()) {
    for (auto c : children_)
      c->removeTris(ixs);
  } else {
    // removing and re-indexing
    // @@ the following reindexing could probably be more efficiently implemented
    vector<uint> new_indices_;
    for (uint i = 0; i != (uint)indices_.size(); ++i) {
      const uint cur_ix = indices_[i];
      const auto f = find_if(ixs.begin(), ixs.end(), [cur_ix](uint u) { return u >= cur_ix;});

      if (f != ixs.end() && cur_ix == *f)
        // index is among those to be removed.  Do not insert it into the new vector
        continue;

      const uint shift = (uint)(f - ixs.begin());
      new_indices_.push_back(cur_ix - shift);
    }
    indices_.swap(new_indices_);
  }
  if (cur_depth_ > 0)
    return;
  
  // we are at top level, so we will remove the actual triangles as well
  vector<uint> tmp(tris_->size(), 1);
  for (auto i : ixs)
    tmp[i] = 0;

  shared_ptr<vector<Triangle>> tris_new(new vector<Triangle>());
  for (uint i = 0; i != (uint)tris_->size(); ++i) 
    if (tmp[i])
      tris_new->push_back((*tris_)[i]);

  tris_->swap(*tris_new);
    
}
  
// ----------------------------------------------------------------------------
void TriangleOctTree::include_last_triangle()
// ---------------------------`-------------------------------------------------
{
  if (is_subdivided()) {
    const auto included_octs = determine_octs(tris_->back());
    for (uint i = 0; i != 8; ++i)
      if (included_octs[i]) 
        children_[i]->include_last_triangle();
  } else {
    indices_.push_back(uint(tris_->size()-1));
    reorganize_if_necessary();
  }
}

// ----------------------------------------------------------------------------  
bool TriangleOctTree::is_subdivided() const
// ----------------------------------------------------------------------------
{
  return (bool) children_[0];
}
  
// ============================================================================
void TriangleOctTree::getIntersectionCandidates(const std::array<Point3D, 3>& tri,
                                                vector<uint>& candidates,
                                                bool clear) const
// ============================================================================
{
  candidates.resize(tris_->size());
  if (clear)
    fill(candidates.begin(), candidates.end(), 0);
  
  if (is_subdivided()) {
    // collect relevant indices from each of the children
    const auto octs = determine_octs(tri);
    for (uint i = 0; i != 8; ++i) 
      if (octs[i])
        children_[i]->getIntersectionCandidates(tri, candidates, false);
  } else {
    for (auto ix : indices_)
      candidates[ix] = 1;
  }
}

// ----------------------------------------------------------------------------
array<bool, 8> TriangleOctTree::determine_octs(const array<Point3D, 3>& t) const
// ----------------------------------------------------------------------------
{
  array<bool, 8> result;
  determine_octs(t, result);
  return result;
}

  // ----------------------------------------------------------------------------
array<bool, 8> TriangleOctTree::determine_octs(const Triangle& t) const
// ----------------------------------------------------------------------------
{
  array<bool, 8> result;
  array<Point3D, 3> tripoints {points_[t[0]], points_[t[1]], points_[t[2]]};
  
  determine_octs(tripoints, result);
  return result;
}

// ----------------------------------------------------------------------------    
void TriangleOctTree::determine_octs(const Triangle& t,
                                     std::array<bool, 8>& result) const
// ----------------------------------------------------------------------------    
{
  array<Point3D, 3> tripoints {points_[t[0]], points_[t[1]], points_[t[2]]};
  determine_octs(tripoints, result);
}
  
// ----------------------------------------------------------------------------  
void TriangleOctTree::determine_octs(const array<Point3D, 3>& t,
                                     array<bool, 8>& result) const
// ----------------------------------------------------------------------------
{
  fill(result.begin(), result.end(), true);

  // eliminate octs that are guaranteed not to contain any part of the triangle
  uint sum[3];
  for (uint dim = 0; dim != 3; ++dim) 
    sum[dim] = (t[0][dim] < midvals_[dim]) +
               (t[1][dim] < midvals_[dim]) +
               (t[2][dim] < midvals_[dim]);

  for (uint i = 0; i != 8; ++i) {
    const uint xi = i%2;     // zero or one
    const uint yi = (i/2)%2; // zero or one
    const uint zi = (i/4);   // zero or one

    if ((sum[0] == (xi ? 3 : 0)) ||
        (sum[1] == (yi ? 3 : 0)) ||
        (sum[2] == (zi ? 3 : 0)))
      result[i] = false;
  }
}
  
// ----------------------------------------------------------------------------
void TriangleOctTree::reorganize_if_necessary()
// ----------------------------------------------------------------------------
{
  if ( ((uint)indices_.size() > max_num_ixs_) && (cur_depth_ < max_depth_ - 1)) {
    //cout << "Number: " << indices_.size() << endl;
    assert(!is_subdivided()); // indices_ should be empty otherwise

    // construct children
    for (uint i = 0; i != 8; ++i) 
      children_[i] = shared_ptr<TriangleOctTree>
        (new TriangleOctTree(points_, tris_, sub_box(bbox_, i), max_num_ixs_,
                             cur_depth_ + 1, max_depth_));
    
    // distribute triangles
    array<bool, 8> octs;
    for (const auto ix : indices_) {
      determine_octs((*tris_)[ix], octs);
      // cout << "Octs: "; for (uint j = 0; j != 8; ++j) { cout << octs[j] << " ";}; cout << endl;
      for (uint i = 0; i != 8; ++i) 
        if (octs[i])
          (children_[i]->indices_).push_back(ix);
    }
    // for (uint i = 0; i != 8; ++i) {
    //   if (children_[i])
    //     cout << "Child " << i << " contains " << children_[i]->indices_.size() << endl;
    // }
    // if (children_[0] && children_[1] && children_[2] && children_[3] &&
    //     children_[4] && children_[5] && children_[6] && children_[7]) {
    //   if ( (children_[0]->indices_.size() == 9) &&
    //        (children_[1]->indices_.size() == 9) &&
    //        (children_[2]->indices_.size() == 25)  &&
    //        (children_[3]->indices_.size() == 17)) {
    //     double krull = 0;
    //   }
      
    // }
    
    indices_.clear();  // will not be used anymore

    // organize children if necessary
    for (uint i = 0; i != 8; ++i)
      if (children_[i]) 
        children_[i]->reorganize_if_necessary();

    //test_integrity();
  }
}

// ----------------------------------------------------------------------------
void TriangleOctTree::test_integrity()
// ----------------------------------------------------------------------------
{
  if (is_subdivided()) {
    for (uint i = 0; i != 8; ++i)
      if (children_[i])
        children_[i]->test_integrity();
  } else {
    const double EPS = numeric_limits<double>::epsilon();
    for (auto ix : indices_) {
      const Triangle& t = (*tris_)[ix];
      array<bool, 3> found {false, false, false};
      array<bool, 3> less  {false, false, false};
      array<bool, 3> more  {false, false, false};
      for (uint i = 0; i != 3; ++i) {
        const Point3D& p = points_[t[i]];
        for (uint dim = 0; dim != 3; ++dim) 
          if ( (p[dim] + EPS >= bbox_[2*dim]) && p[dim] - EPS <= bbox_[2*dim+1])
            found[dim] = true;
          else if (p[dim] < bbox_[2*dim])
            less[dim] = true;
          else if (p[dim] > bbox_[2*dim+1])
            more[dim] = true;
      }
      for (uint dim = 0; dim != 3; ++dim)
        found[dim] = found[dim] || (more[dim] && less[dim]);
                           
      assert(found[0] && found[1] && found[2]);
    }
  }
}

// ----------------------------------------------------------------------------
TriangleOctTree::TriangleOctTree(const Point3D* const pts,
                                 shared_ptr<vector<Triangle>> tris,
                                 const std::array<double, 6>& bbox,
                                 const uint max_num_ixs,
                                 const uint cur_depth,
                                 const uint max_depth)
// ----------------------------------------------------------------------------
  : max_num_ixs_(max_num_ixs),
    cur_depth_(cur_depth),
    max_depth_(max_depth),
    points_(pts),
    tris_(tris),
    bbox_(bbox),
    midvals_(compute_midvals(bbox_))
{}

 
  
}; // end namespace TesselateUtils
