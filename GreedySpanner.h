//
// Created by FNU Shariful on 6/7/22.
//
#include <iostream>
#include "vector"
#include <CGAL/Simple_cartesian.h>
#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "CGALComponents.h"
#include "DelaunayTriangulation.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
//typedef pair<unsigned long,unsigned long > Edge;

typedef size_t index_t;
typedef size_t cone_t;
//typedef variant<index_t,number_t> mixed_t;

typedef pair<index_t, index_t> index_tPair;
typedef index_tPair Edge;

#ifndef TEST_GREEDYSPANNER_H
#define TEST_GREEDYSPANNER_H

//using namespace std;

// FG-greedy with Delaunay; much faster than orginal greedy
void constructFG_GreedySpannerV4(vector<Point> &P, vector<Edge> &spannerEdges, const double t, unordered_map<int, vector<int>> &G) {

  vector<vector<double>> distances(P.size(), vector<double>(P.size(), 0)),
      shortestPathLength(P.size(), vector<double>(P.size(), DBL_MAX));
  spannerEdges.reserve(3 * P.size() );
//  vector<unordered_set<unsigned>> G(P.size());

  // Using Del triang. as the starting point is somewhat helpful; slightly faster but puts a bit more edges
  unordered_map<Point,unsigned> idMap;
  unsigned index = 0;
  for(Point p : P)
    idMap[p] = index++;

//  vector<Edge> delaunayEdges;
  DelaunayTriangulation DT(P,spannerEdges,idMap);
//
//  for( Edge e : delaunayEdges ) {
//    G[idMap[e.first]].insert(idMap[e.second]);
//    G[idMap[e.second]].insert(idMap[e.first]);
//    spannerEdges.emplace_back(e.first, e.second);
//  }
  /////////////////// comment: code runs much faster if sqrt in distance calculation is not used but gives wrong s.f.

  vector<Edge> completeGraphE;
  completeGraphE.reserve((P.size() * (P.size()-1))/2);
  for (unsigned i = 0; i < P.size() - 1; i++)
    for (unsigned j = i + 1; j < P.size(); j++) {
      completeGraphE.emplace_back(make_pair(i, j));
      distances[i][j] = distances[j][i] = std::sqrt(squared_distance(P[i], P[j])) ;
    }

  sort(completeGraphE.begin(), completeGraphE.end(), [distances](const Edge &e1, const Edge &e2) {
    return distances[e1.first][e1.second] < distances[e2.first][e2.second];
  });

  /////////// boost's stuff
  using namespace boost;
  typedef adjacency_list<vecS, vecS, undirectedS, no_property,
                         property<edge_weight_t, double, property<edge_weight2_t, double> > > Graph;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  Graph g(&spannerEdges[0], &spannerEdges[spannerEdges.size()], P.size());
  vector<double> d(num_vertices(g));
  vector<vertex_descriptor> p(num_vertices(g));
  property_map<Graph, edge_weight_t>::type w = get(edge_weight, g);
  vector<double> weightsForDijkstraGraph(spannerEdges.size());
  weightsForDijkstraGraph.reserve(6 * P.size());
  graph_traits<Graph>::edge_iterator edge, e_end;
  boost::tie(edge, e_end) = edges(g);

  for (unsigned i = 0; i < spannerEdges.size(); i++)
    weightsForDijkstraGraph[i] = distances[spannerEdges[i].first][spannerEdges[i].second];

  double *wp = &weightsForDijkstraGraph[0];
  //////////////////////

  for (const Edge &e: completeGraphE) {
//    if (G[e.first].find(e.second) != G[e.first].end()) // if the edge is present, do nothing
//      continue;

    if(find(G[e.first].begin(), G[e.first].end(),e.second) != G[e.first].end() )
      continue ;

    double tTimesdistance = t * distances [e.first] [e.second] ;

    if (shortestPathLength[e.first][e.second] <= tTimesdistance)
      continue;
    /////////////////// using boost's dijkstra
    for (; edge != e_end; ++edge)
      w[*edge] = *wp++;

    dijkstra_shortest_paths(g, vertex(e.first, g), boost::predecessor_map(&p[0]).distance_map(&d[0]));

    for (unsigned i = 0; i < d.size(); i++)
      shortestPathLength[e.first][i] = shortestPathLength[i][e.first] = std::min(shortestPathLength[e.first][i], d[i]);
    ////////////////////
    if (shortestPathLength[e.first][e.second] > tTimesdistance ) {
      spannerEdges.emplace_back(e.first, e.second);
      G[e.first].emplace_back(e.second);
      G[e.second].emplace_back(e.first);
      boost::add_edge(e.first, e.second, distances[e.first][e.second], g);
      weightsForDijkstraGraph.emplace_back(distances[e.first][e.second]);
    }
  }
//  return spannerEdges;
}


#endif // TEST_GREEDYSPANNER_H
