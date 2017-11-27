#include <map>
#include <list>
#include <set>
#include <algorithm>
#include "tritools.h"

#include "ttl/ttl.h"
#include "ttl/halfedge/HeDart.h"
#include "ttl/halfedge/HeTraits.h"
#include "ttl/halfedge/HeTriang.h"



using namespace Go;
using namespace std;
using namespace ttl;
using namespace hed;
using namespace TriTools;

namespace {

// ----------------------------------------------------------------------------
inline double cross_value (const Point& p1, const Point& p2)
// ----------------------------------------------------------------------------  
{
  return p1[0]*p2[1] - p2[0] * p1[1];
}
  
// ----------------------------------------------------------------------------
// this function assumes that the line equation is on normalized form, i.e. a^2 + b^2 = 1,
// where a = line_equation[0] and b = line_equation[1].
inline double signed_distance(const array<double, 3>& line_equation, const Point& p)
// ----------------------------------------------------------------------------  
{
  return line_equation[0] * p[0] + line_equation[1] * p[1] + line_equation[2];
}

// ----------------------------------------------------------------------------  
inline array<double, 4> get_tribox(const Point& p1, const Point& p2, const Point& p3)
// ----------------------------------------------------------------------------  
{
  return
    { min(min(p1[0], p2[0]), p3[0]),  // xmin
      max(max(p1[0], p2[0]), p3[0]),  // xmax
      min(min(p1[1], p2[1]), p3[1]),  // ymin 
      max(max(p1[1], p2[1]), p3[1])   // ymax
    }; 
}

// ----------------------------------------------------------------------------  
Triangulation create_delaunay_from_loops(const vector<vector<Point>>& loops)
// ----------------------------------------------------------------------------
{
  vector<Node*> nodes;
  for (const auto vec : loops) // loop over loops
    for (const auto pt : vec) // loop over points in loop
      nodes.emplace_back(new Node {pt[0], pt[1]});
  sort(nodes.begin(), nodes.end(), [](const Node* n1, const Node* n2) {
      return (n1->x() < n2->x() || (n1->x() == n2->x() && n1->y() < n2->y()));
    });
  Triangulation result;
  result.createDelaunay(nodes.begin(), nodes.end());
  return result;
}

// ----------------------------------------------------------------------------
Edge* get_leading_edge(Edge* e) {
// ----------------------------------------------------------------------------
  if (!e) {
    return NULL;
  } else if (e->isLeadingEdge()) {
    return e;
  } else if (e->getNextEdgeInFace()->isLeadingEdge()) {
    return e->getNextEdgeInFace();
  }
  assert(e->getNextEdgeInFace()->getNextEdgeInFace()->isLeadingEdge());
  return e->getNextEdgeInFace()->getNextEdgeInFace();
}
  
// ----------------------------------------------------------------------------
void recursively_find_triangles(Edge* const e, set<Edge*>& result)
// ----------------------------------------------------------------------------
{
  result.insert(get_leading_edge(e));
  
  array<Edge*, 3> neighs;
  neighs[0] = get_leading_edge(e->getTwinEdge());
  neighs[1] = get_leading_edge(e->getNextEdgeInFace()->getTwinEdge());
  neighs[2] = get_leading_edge(e->getNextEdgeInFace()->getNextEdgeInFace()->getTwinEdge());

  for (auto n : neighs)
    if (bool(n) & (result.find(n) == result.end()))
	recursively_find_triangles(n, result);
}
  
// ----------------------------------------------------------------------------
template<typename PointIterator>
void impose_boundary(Triangulation& tri,
		     PointIterator pts_begin,
		     PointIterator pts_end)
// ----------------------------------------------------------------------------
{
  // mapping coordinates of existing nodes to representative edges in triangulation
  map<Point, Edge*> point_to_edge;
  const list<Edge*>* const edges = tri.getEdges();
  for (const auto e : *edges) {
    const Node* const n = e->getSourceNode();
    point_to_edge[Point(n->x(), n->y())] = e;
  }

  // insert missing points
  for (auto pt = pts_begin; pt != pts_end; ++pt)
    if (point_to_edge.find(*pt) == point_to_edge.end()) {
      // the point is not already present, so we should insert it
      auto d = tri.createDart();
      insertNode<TTLtraits>(d, *new Node((*pt)[0], (*pt)[1]));
      point_to_edge[*pt] = d.getEdge();
    }

  // set constraints here
  vector<Edge*> outside_edges;
  for (auto cur = pts_begin; cur != pts_end; ++cur) {
    Edge* e1 = point_to_edge[*cur];
    Edge* e2 = point_to_edge[ (cur+1 != pts_end) ? *(cur+1) : *pts_begin];
    Dart d1(e1), d2(e2);
    auto dart = insertConstraint<TTLtraits>(d1, d2, true);
    dart.getEdge()->setConstrained();

    // if the constriant is an interior edge, make it a boundary edge
    if (!isBoundaryEdge(dart)) 
      outside_edges.push_back(dart.getEdge());
  }

  // setting interior constrained edges to boundary edges
  for (auto e : outside_edges)
    if (e) {
      auto twin_edge = e->getTwinEdge();
      e->setTwinEdge(NULL);
      twin_edge->setTwinEdge(NULL);
    }
  
  // identify all triangles outside the introduced boundary
  set<Edge*> triangles_found;
  for (auto e : outside_edges)
    if (e)
      recursively_find_triangles(e, triangles_found);

  // delete triangles
  for (auto t : triangles_found)
    tri.removeTriangle(*t);
    
}

// ----------------------------------------------------------------------------  
vector<pair<array<Point, 3>, array<bool, 3>>>
export_triangles(const Triangulation& tri)
// ----------------------------------------------------------------------------
{
  vector<pair<array<Point, 3>, array<bool, 3>>> result;
  for (auto edge : tri.getLeadingEdges()) {
    const auto e1 = edge;
    const auto e2 = e1->getNextEdgeInFace();
    const auto e3 = e2->getNextEdgeInFace();
    const Node* const n1 = e1->getSourceNode();
    const Node* const n2 = e2->getSourceNode();
    const Node* const n3 = e3->getSourceNode();

    result.emplace_back(pair<array<Point, 3>, array<bool, 3>>
			{{Point(n1->x(), n1->y()),
	                  Point(n2->x(), n2->y()),
	                  Point(n3->x(), n3->y())},
	                 {e1->isConstrained(),
			  e2->isConstrained(),
			  e3->isConstrained()}});
  }
  return result;
}

// ----------------------------------------------------------------------------
inline bool inside_box(const array<double, 4>& box, const Point& p)
// ----------------------------------------------------------------------------
{
  return (p[0] > box[0] && p[0] < box[1] &&
	  p[1] > box[2] && p[1] < box[3]);
}

// ----------------------------------------------------------------------------
inline bool segment_close(const Point& p1, const Point& p2,
			  const Point& p, double margin)
// ----------------------------------------------------------------------------
{
  // testing against endpoints
  if (min(p1.dist2(p), p2.dist2(p)) < margin*margin)
    return true;

  // testing against segment
  return ( ( (p2-p1) * (p-p1) > 0 ) &&
	   ( (p1-p2) * (p-p2) > 0 ) &&
	   ( abs(signed_distance(TriTools::normalized_line_equation(p1,p2), p)) < margin));
}

// ----------------------------------------------------------------------------
void remove_points_close_to_segment(const Point& p1, const Point& p2,
				    const double margin,
				    const vector<Point>& points,
				    vector<int>& result)
// ----------------------------------------------------------------------------
{
  // compute bounding box
  array<double, 4> bbox {min(p1[0], p2[0])-margin, // xmin
                         max(p1[0], p2[0])+margin, // xmax
                         min(p1[1], p2[1])-margin, // xmin
                         max(p1[1], p2[1])+margin}; // xmax
  for (size_t i = 0; i != points.size(); ++i) 
    if (result[i] && inside_box(bbox, points[i]))
      if (segment_close(p1, p2, points[i], margin))
	result[i] = 0;
}
  
// ----------------------------------------------------------------------------
void remove_points_close_to_boundaries(const vector<vector<Point>>& loops,
				       const double margin,
				       const vector<Point>& points,
				       vector<int>& result)
// ----------------------------------------------------------------------------
{
  for (auto l : loops) 
    for (size_t pt_ix = 0; pt_ix != l.size(); ++pt_ix) {
      const Point& p1 = l[pt_ix];
      const Point& p2 = l[((pt_ix+1) != l.size() ? pt_ix+1 : 0)];
      remove_points_close_to_segment(p1, p2, margin, points, result);
    }
}

// ----------------------------------------------------------------------------
pair<vector<Node*>, vector<Node*>>
get_interior_and_boundary_nodes(const Triangulation& tri)
// ----------------------------------------------------------------------------
{
  const list<Edge*> all_edges = *(tri.getEdges());
  vector<Node*> boundary_nodes, interior_nodes;
  for (auto e : all_edges)
    if (ttl::isBoundaryEdge(Dart(e)))
      boundary_nodes.push_back(e->getSourceNode());
    else
      interior_nodes.push_back(e->getSourceNode());

  sort(interior_nodes.begin(), interior_nodes.end());
  auto last = unique(interior_nodes.begin(), interior_nodes.end());
  interior_nodes.erase(last, interior_nodes.end());

  return {interior_nodes, boundary_nodes};
}

// ----------------------------------------------------------------------------
SurfaceTriangulation convert_to_surface_triangulation(const Triangulation& tri)
// ----------------------------------------------------------------------------
{
  // making mapping from nodes to indices, and converting nodes to vector of
  // Go::Points.
  const auto int_and_bnd_nodes = get_interior_and_boundary_nodes(tri);
  
  map<const Node*, int> node_to_ix;
  vector<Point> points;
  int counter = 0;
  // indexing interior and boundary nodes
  const vector<Node*>& int_nodes = int_and_bnd_nodes.first;
  const vector<Node*>& bnd_nodes = int_and_bnd_nodes.second;
  for (auto it : int_nodes) {
    points.push_back({ it->x(), it->y()});
    node_to_ix[it] = counter++;
  }
  for (auto it : bnd_nodes) {
    points.push_back({ it->x(), it->y()});
    node_to_ix[it] = counter++;
  }
  
  // constructing topology
  vector<array<int, 3>> triangles;
  const list<Edge*>& leading_edges = tri.getLeadingEdges();
  for (auto e = leading_edges.begin(); e != leading_edges.end(); ++e) {
    const Node* n1 = (*e)->getSourceNode();
    const Node* n2 = (*e)->getNextEdgeInFace()->getSourceNode();
    const Node* n3 = (*e)->getNextEdgeInFace()->getNextEdgeInFace()->getSourceNode();
    triangles.push_back({node_to_ix[n1], node_to_ix[n2], node_to_ix[n3]});
  }
  
  // returning result
  return {points, triangles, int_nodes.size()};
}

  
}; // end anonymous namespace 

namespace TriTools {

// ----------------------------------------------------------------------------
array<double, 3> normalized_line_equation(const Point& pt_first,
					  const Point& pt_second)
// ----------------------------------------------------------------------------
{
  const Point dp = pt_second - pt_first;
  const double inv_length = 1.0/dp.length();

  return {
     dp[1] * inv_length,
    -dp[0] * inv_length,
    -cross_value(pt_first, pt_second) * inv_length
  };
}

// ----------------------------------------------------------------------------
vector<int> points_inside_triangle(const Point& p1,
				   const Point& p2,
				   const Point& p3,
				   const vector<Point>& points,
				   const array<double, 3>& eps)
// ----------------------------------------------------------------------------
{
  vector<int> result(points.size(), 0);
  transform(points.begin(), points.end(), result.begin(), [&] (const Point& p) {
      return (int)(
	(signed_distance(normalized_line_equation(p1, p2), p) <= -eps[0]) &
	(signed_distance(normalized_line_equation(p2, p3), p) <= -eps[1]) &
	(signed_distance(normalized_line_equation(p3, p1), p) <= -eps[2]));});
  return result;
}

// ----------------------------------------------------------------------------
vector<int> points_inside_loops(const vector<vector<Point>>& loops,
				const vector<Point>& points,
				const double margin)
// ----------------------------------------------------------------------------
{
  const vector<pair<array<Point, 3>, array<bool, 3>>> triangles =
	       triangulate_boundary(loops);

  vector<int> result(points.size(), 0); // we initially consider all points outside
  for (const auto tri : triangles) {
    auto cur_inside = points_inside_triangle(tri.first[0],
					     tri.first[1],
					     tri.first[2],
					     points, {0l, 0l, 0l});
    transform(result.begin(), result.end(), cur_inside.begin(), result.begin(),
	      [](int a, int b) {return int(a || b);});
  }

  // removing points that are too close to the boundary
  if (margin > 0)
    remove_points_close_to_boundaries(loops, margin, points, result);

  return result;
}

// ----------------------------------------------------------------------------
vector<pair<array<Point, 3>, array<bool, 3>>>
triangulate_boundary(const vector<vector<Point>>& loops)
// ----------------------------------------------------------------------------
{
  // Create delauney triangulation containing all boundary points
  Triangulation triang = create_delaunay_from_loops(loops);

  // auto nodes = *triang.getNodes();
  // for (n : nodes)
  //   cout << "Node is: " << n->x() << " " << n->y() << '\n';
  
  // constrain boundary edges
  for (auto l : loops)
    impose_boundary(triang, l.begin(), l.end());

  // generate triangles of the result
  return export_triangles(triang);
  
}

// ----------------------------------------------------------------------------
// @@ This function may be made more efficient
SurfaceTriangulation
triangulate_with_boundaries(const vector<Point>& points,
			    const vector<vector<Point>>& loops)
// ----------------------------------------------------------------------------
{
  // Create initial delauney triangulation from boundary points
  Triangulation triang = create_delaunay_from_loops(loops);

  // Constrain boundary edges
  for (auto l : loops)
    impose_boundary(triang, l.begin(), l.end());

  // Add the rest of the points (assuming they are all inside the valid region)
  Dart d = triang.createDart();
  for (auto p : points) {
    Node* nd = new Node(p[0], p[1]);
    insertNode<TTLtraits>(d, *nd);
  }
  return convert_to_surface_triangulation(triang);
}

  
}; // end namespace TriTools
