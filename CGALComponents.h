//
// Created by FNU Shariful on 6/7/22.
//

#ifndef TEST_CGALCOMPONENTS_H
#define TEST_CGALCOMPONENTS_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <unordered_set>
#include <vector>

#include <boost/heap/fibonacci_heap.hpp>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/ch_jarvis.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>

#include <QApplication>
#include <QLabel>
#include <QString>
#include <QTranslator>
#include <QtWidgets>

using namespace std;
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Segment_2 Segment;
typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
typedef CGAL::Alpha_shape_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Alpha_shape_2<CGAL::Delaunay_triangulation_2<K, Tds>> Alpha_shape;

typedef CGAL::Search_traits_2<K> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;

typedef pair<unsigned long, unsigned long> Edge;

// class Point : public K::Point_2 {
// public:
//   unsigned long id = -99; // -99 is a dummy value for id
//   string color = "black";
//   long int v = 0, h = 0;
//
//   Point() = default;
//
//   Point(double X, double Y) : K::Point_2(X, Y) { }
//   Point(double X, double Y, unsigned long k) : K::Point_2(X, Y) { id = k; }
//
//   friend ostream &operator<<(ostream &strm, const Point &p) {
//     return strm << p.id << ": (" << p.x() << ", " << p.y() << ")";
//   }
// };

inline double L2distance(const Point &p1, const Point &p2) {
  return std::sqrt(squared_distance(p1, p2));
}

bool turnOnPdfs = false;

void generatePointsInsideASquare(const unsigned n, const double sizeOfSquare,
                                 vector<Point> &P) {
  vector<K::Point_2> points;
  typedef Random_points_in_square_2<K::Point_2,
                                    Creator_uniform_2<double, K::Point_2>>
      Point_generator;
  Point_generator g(sizeOfSquare / 2);

  // random_convex_set_2(n, back_inserter(points), g);

  std::copy_n(g, n, back_inserter(P));

  //  unsigned id = 0;
  //  for( K::Point_2 p : points )
  //    P.emplace_back( Point(p.x(), p.y(), id++) );

  // perturb_points_2( P.begin(), P.end(), 0.0001, 0.0001 );
}

void generatePointsInsideADisc(const unsigned n, const double radiusOfDisk,
                               vector<Point> &P) {
  vector<K::Point_2> points;
  typedef Random_points_in_disc_2<K::Point_2,
                                  Creator_uniform_2<double, K::Point_2>>
      Point_generator;
  Point_generator g(radiusOfDisk);

  std::copy_n(g, n, back_inserter(points));
  // perturb_points_2( points.begin(), points.end(), 0.0001, 0.0001 );

  unsigned id = 0;
  for (K::Point_2 p : points)
    P.emplace_back(Point(p.x(), p.y(), id++));
}

void generatePointsInsideASquareNormal(const unsigned pointsInACuster,
                                       const unsigned numberOfClusters,
                                       vector<Point> &P) {
  vector<K::Point_2> points;

  std::mt19937 rngX(std::random_device{}());
  std::mt19937 rngY(std::random_device{}());
  std::default_random_engine generatorX(rngX()), generatorY(rngY());
  std::normal_distribution<double> distributionX(2.0, 2.0),
      distributionY(2.0, 2.0);

  // std::poisson_distribution<int> distributionX(4.1), distributionY(4.1);

  std::mt19937 rngShift(std::random_device{}());

  std::uniform_int_distribution shiftDistribution(0, INT32_MAX);

  unsigned shiftX, shiftY;
  unordered_set<pair<unsigned, unsigned>, boost::hash<pair<unsigned, unsigned>>>
      S;

  for (unsigned c = 0; c < numberOfClusters; c++) {
    if (c != 0) {
      shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
      shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);

      while (!(S.find(make_pair(shiftX, shiftY)) == S.end())) {
        shiftX = shiftDistribution(rngShift) % (20 * numberOfClusters);
        shiftY = shiftDistribution(rngShift) % (20 * numberOfClusters);
      }
    } else
      shiftX = shiftY = 0;

    S.insert(make_pair(shiftX, shiftY));

    for (unsigned i = 0; i < pointsInACuster; i++) {
      double x = distributionX(generatorX) + shiftX;
      double y = distributionY(generatorY) + shiftY;
      points.emplace_back(K::Point_2(x, y));
    }
  }

  // perturb_points_2( points.begin(), points.end(), 0.0001, 0.0001 );

  unsigned id = 0;
  for (Point p : points)
    P.emplace_back(Point(p.x(), p.y(), id++));

  //    for(const Point &p : P)
  //      printf("%ld: (%.10lf, %.10lf)\n",p.id,p.x(),p.y());
}

void generateContiguousPointsOnAGrid(const unsigned n, vector<Point> &P) {
  vector<Point> points;
  points_on_square_grid_2(10, n, std::back_inserter(points),
                          Creator_uniform_2<int, Point>());
  // perturb_points_2( points.begin(), points.end(), 0.0001, 0.0001 );

  // sort(points.begin(),points.end());

  unsigned id = 0;
  for (Point p : points)
    P.emplace_back(Point(p.x(), p.y(), id++));
  //    for(const Point &p : P)
  //        printf("%ld: (%.10lf, %.10lf)\n",p.id,p.x(),p.y());
}

void generateRandomPointsOnAGrid(const int n, vector<Point> &P) {

  vector<K::Point_2> points;

  unordered_set<pair<int, int>, boost::hash<pair<int, int>>> S;
  std::mt19937 rngX(std::random_device{}());
  std::mt19937 rngY(std::random_device{}());
  std::uniform_int_distribution xDistribution(0, (int)ceil(0.7 * n)),
      yDistribution(0, (int)ceil(0.7 * n));

  int count = 0;

  while (count < n) {
    int x = xDistribution(rngX), y = yDistribution(rngY);

    if (S.find(make_pair(x, y)) == S.end()) {
      points.emplace_back(Point(x, y, count));
      S.insert(make_pair(x, y));
      count++;
    }
  }

  // perturb_points_2(points.begin(), points.end(), 0.0001, 0.0001 );

  unsigned id = 0;
  for (K::Point_2 p : points)
    P.emplace_back(Point(p.x(), p.y(), id++));

  //    sort(points.begin(),points.end());
  //    for(const Point &p : P)
  //        printf("%ld: (%.10lf, %.10lf)\n",p.id,p.x(),p.y());
}

void generateRandomInsideAnAnnulus(const unsigned n, const double r2,
                                   const double r1, vector<Point> &P) {

  assert(r2 > r1);
  assert(n > 1);

  vector<K::Point_2> points;

  // std::default_random_engine generator;
  std::mt19937 generator(std::random_device{}());

  std::uniform_real_distribution<double> distributionR(r1, r2),
      distributionT(0, 1);

  for (unsigned i = 0; i < n; i++) {
    double t = 2 * M_PI * distributionT(generator);
    double r = distributionR(generator);
    points.emplace_back(K::Point_2(r * cos(t), r * sin(t)));
  }

  // perturb_points_2(points.begin(), points.end(), 0.0001, 0.0001 );

  unsigned id = 0;
  for (K::Point_2 p : points)
    P.emplace_back(Point(p.x(), p.y(), id++));
}

void drawSpannerUsingQT(const std::vector<Point> &Vertices,
                        const vector<Edge> &egdes,
                        const bool labels) {
  //    assert( !polygonVertices.empty() );

  const double pointSize = 2; // adjust this value to change the size of the points
  /***************************************************/
  QPicture pi;
  QPainter canvas(&pi);
  canvas.setRenderHint(QPainter::Antialiasing);
  // canvas.setFont(QFont("Times", 12));
  //  DO NOT TOUCH THESE Qt STATEMENTS!!!
  /***************************************************/

  canvas.setBrush(Qt::white);

  std::vector<QPointF> verticesForQTPolygon;
  verticesForQTPolygon.reserve(Vertices.size());
  for (Point p : Vertices)
    verticesForQTPolygon.emplace_back(QPointF(p.x(), p.y()));

  canvas.drawPolygon(&verticesForQTPolygon[0],
                     (int)verticesForQTPolygon.size());

  for (Edge e : egdes) {
    QPointF endPointA(Vertices[e.first].x(),
                      Vertices[e.first].y()),
        endPointB(Vertices[e.second].x(), Vertices[e.second].y());
    canvas.drawLine(endPointA, endPointB);
  }

  unsigned id = 0;
  for (Point p : Vertices) {

//    if (vertexColors[id] == Qt::red) {
//      canvas.setBrush(Qt::red);
//      canvas.setPen(Qt::red);
//    } else if (vertexColors[id] == Qt::darkGreen) {
//      canvas.setBrush(Qt::darkGreen);
//      canvas.setPen(Qt::darkGreen);
//    } else {
//      canvas.setBrush(Qt::blue);
//      canvas.setPen(Qt::blue);
//    }

    canvas.drawEllipse(QPointF(p.x(), p.y()), pointSize, pointSize);
    if (labels) {
      canvas.setBrush(Qt::black);
      canvas.setPen(Qt::black);
      QPointF textPos(p.x() + 4.0, p.y() + 4.0);
      canvas.drawText(textPos, QString(std::to_string(id).c_str()));
    }
    id++;
  }

  /*********************************************/
  canvas.end();
  auto l = new QLabel;
  l->setStyleSheet("QLabel { background-color : white; color : black; }");
  l->setPicture(pi);
  l->setWindowTitle("Polygon Triangulation");
  l->show();
  // l->showMaximized();
  QApplication::exec();
  // DO NOT TOUCH THESE Qt STATEMENTS!!!
  /***************************************************/
}

#endif // TEST_CGALCOMPONENTS_H
