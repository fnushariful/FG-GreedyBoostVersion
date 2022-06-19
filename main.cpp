#include <iostream>
#include <CGAL/Simple_cartesian.h>
#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//#include <CGAL/Timer.h>
#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

#include "GreedySpanner.h"
#include "CGALComponents.h"
#include "DelaunayTriangulation.h"

//#include <QtWidgets>
//#include <QApplication>
//#include <QLabel>
//#include <QString>
//#include <QTranslator>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    for( int i=1; i<=100; i++ ) {
      vector<Point> P;
      unsigned numberOfPoints = i*1000;
      unsigned sizeOfSquare = 500;
      generatePointsInsideASquare(numberOfPoints, sizeOfSquare, P);

      double stretchFactor = 1.1;
      unordered_map<int, vector<int>> adjMap, adjMap2;
      vector<Edge> spannerEdges;

      CGAL::Real_timer clock;
      clock.start();
      constructGreedySpannerV5(P, spannerEdges, stretchFactor, adjMap);
      drawSpannerUsingQT(P,spannerEdges,false);
      cout << spannerEdges.size() << endl;
      clock.stop();
      cout << "Time Needed " << clock.time() << " seconds = " << clock.time() / 60
           << " minutes" << endl;
    }



//    Point_2 p(1,1), q(10,10);
//    std::cout << "p = " << p << std::endl;
//    std::cout << "q = " << q.x() << " " << q.y() << std::endl;
//    std::cout << "sqdist(p,q) = "
//              << CGAL::squared_distance(p,q) << std::endl;
//    Segment_2 s(p,q);
//    Point_2 m(5, 9);
//    std::cout << "m = " << m << std::endl;
//    std::cout << "sqdist(Segment_2(p,q), m) = "
//              << CGAL::squared_distance(s,m) << std::endl;
//    std::cout << "p, q, and m ";
//    switch (CGAL::orientation(p,q,m)){
//        case CGAL::COLLINEAR:
//            std::cout << "are collinear\n";
//            break;
//        case CGAL::LEFT_TURN:
//            std::cout << "make a left turn\n";
//            break;
//        case CGAL::RIGHT_TURN:
//            std::cout << "make a right turn\n";
//            break;
//    }
//    std::cout << " midpoint(p,q) = " << CGAL::midpoint(p,q) << std::endl;

    return 0;
}