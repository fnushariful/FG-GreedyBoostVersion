//
// Created by Ghosh, Anirban on 5/6/21.
//

#ifndef NAIVEWSPD_DELAUNAYTRIANGULATION_H
#define NAIVEWSPD_DELAUNAYTRIANGULATION_H

#include "CGALComponents.h"

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "CGAL/spatial_sort.h"

class DelaunayTriangulation {
    typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;

    Delaunay T;
public:
    DelaunayTriangulation( vector<Point> &P, vector<Edge> &E, unordered_map<Point,unsigned > idMap) {

        CGAL::spatial_sort(P.begin(),P.end());

        vector<pair<Delaunay::Point, unsigned> > points;
//
        unsigned index = 0;
        for( Point p : P )
            points.emplace_back(make_pair(Delaunay::Point(p.x(),p.y()),index++));

        T.insert(points.begin(), points.end());

        for(Delaunay::Finite_edges_iterator it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it) {
            Delaunay::Edge e=*it;
            unsigned i1 = e.first->vertex( (e.second+1)%3 )->info();
            unsigned i2 = e.first->vertex( (e.second+2)%3 )->info();
            assert(i1 < P.size() && i2 < P.size());
            E.emplace_back(make_pair(i1,i2));
        }
    }

    DelaunayTriangulation( vector<Point> &P, unordered_map<Point,unsigned > idMap) {
        vector<pair<Delaunay::Point, unsigned> > points;

        for( Point p : P )
            points.emplace_back(make_pair(Delaunay::Point(p.x(),p.y()),idMap[p]));

        T.insert(points.begin(), points.end());
    }

    unsigned findClosestPoint(const Point &p) {
        auto handle = T.nearest_vertex(p);
        return handle->info();
    }
};
#endif //NAIVEWSPD_DELAUNAYTRIANGULATION_H
