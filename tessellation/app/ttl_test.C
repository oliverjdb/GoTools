#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <assert.h>

#include "ttl/ttl.h"
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>

using namespace std;
using namespace hed;
using namespace ttl;

namespace {

  Edge* get_leading_edge(Edge* e) {
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

  void recursively_find_triangles(Edge* const e, set<Edge*>& result)
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

  void impose_boundary(Triangulation& tri,
		       const double* const xy_begin,
		       const double* const xy_end)
  {
    // mapping coordinates of existing nodes to representative edges in
    // triangulation
    typedef pair<double, double> Coord2D;
    map<Coord2D, Edge*> node_to_edge;
    const list<Edge*>* const edges = tri.getEdges();
    for (const auto e : *edges) {
      Node* n = e->getSourceNode();
      node_to_edge[{n->x(), n->y()}] = e;
    }
  
    // making list of coords constituting the boundary
    vector<Coord2D> all_coords;
    assert((xy_end - xy_begin) % 2 == 0);
    for (auto it = xy_begin; it != xy_end; it += 2)
      all_coords.push_back({*it, *(it+1)});
  
    // insert missing nodes
    for (auto c : all_coords) {
      auto map_iter = node_to_edge.find(c);
      if (map_iter == node_to_edge.end()) {
	// not found, so we should insert it
	auto d = tri.createDart();
	insertNode<TTLtraits>(d, *new Node(c.first, c.second));
	node_to_edge[c] = d.getEdge();
      }
    }
    // set constraints here
    vector<Edge*> outside_edges; // former twins
    for (size_t i = 0; i != all_coords.size(); ++i) {
      Edge* e1 = node_to_edge[all_coords[i]];
      Edge* e2 = node_to_edge[all_coords[(i+1) % all_coords.size()]];
      Dart d1(e1), d2(e2);
      auto dart = insertConstraint<TTLtraits>(d1, d2, true);
      dart.getEdge()->setConstrained();
      
      // if the constraint is an interior edge, make it a boundary edge
      if (!isBoundaryEdge(dart)) {
	auto twin_edge = dart.getEdge()->getTwinEdge();
	dart.getEdge()->setTwinEdge(NULL);
	twin_edge->setTwinEdge(NULL);
	outside_edges.push_back(dart.getEdge());
      }
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
};
    
int main()
{
  // Data (node coordinates)
  const vector<double> points {0,0,  4,0,  7,0,
                               7,2,  3,2,  3,6,
                               5,6,  7,6,  7,8,
                               3,8,  0,8,  0,5,  0,2};

  // defining and lexicographically sorting nodes
  vector<Node*> nodes;
  for (size_t i = 0; i != points.size(); i += 2)
    nodes.emplace_back(new Node {points[i], points[i+1]});

  vector<Node*> unsorted_nodes = nodes;
  sort(nodes.begin(), nodes.end(), [](const Node* p1, const Node* p2) {
      return (p1->x() < p2->x()) || (p1->x() == p2->x() && p1->y() < p2->y());
    });

  // Make the triangulation
  Triangulation triang;
  triang.createDelaunay(nodes.begin(), nodes.end());

  //constrain boundary edges
  impose_boundary(triang, &points[0], &points[0] + points.size());

  // const list<Edge*>* edges = triang.getEdges();
  // for (size_t i = 0; i != unsorted_nodes.size(); ++i) {
  //   auto e1 = *find_if(edges->begin(), edges->end(), [&](Edge* e) {
  // 	auto n1 = e->getSourceNode();
  // 	auto n2 = unsorted_nodes[i];
  // 	return (n1->x() == n2->x()) && (n1->y() == n2->y());});
  //   auto e2 = *find_if(edges->begin(), edges->end(), [&](Edge* e) {
  // 	auto n1 = e->getSourceNode();
  // 	auto n2 = unsorted_nodes[(i+1) % unsorted_nodes.size()];
  // 	return (n1->x() == n2->x()) && (n1->y() == n2->y());});

  //   Dart d1(e1), d2(e2);
  //   auto dart = insertConstraint<TTLtraits>(d1, d2, true);
  //   // insertNode<TTLtraits>(d2, *unsorted_nodes[i+1]);
  //   // auto dart = insertConstraint<TTLtraits>(d1, d2, true);
  //   dart.getEdge()->setConstrained();
  // }

  // removing triangles outside constrained region
  
  // while (true) {
  //   bool removed = false;
  //   list<Edge*>* elist = triang.getEdges();
  //   for (auto e : *elist) {
  //     if (isBoundaryEdge(Dart(e)) && !e->isConstrained()) {
  // 	triang.removeTriangle(*e);
  // 	removed = true;
  //     }
  //   }
  //   if (!removed)
  //     break;
  // }
  
  
  // insertNode<TTLtraits>(d1, n1);
  // insertNode<TTLtraits>(d2, n2);
  // auto dart = insertConstraint<TTLtraits>(d1, d2, true);
  // dart.getEdge()->setConstrained();
  
  // plotting result
  ofstream os("triang.dat");
  triang.printEdges(os);
  os.close();
  
  return 0;
};
