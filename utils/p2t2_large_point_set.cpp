#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT> PDT;

typedef PDT::Point          Point;
typedef PDT::Face_handle                                    Face_handle;
typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Face                                           Face;
typedef PDT::Finite_faces_iterator                          Face_iterator;
typedef PDT::Vertex_iterator                               Vertex_iterator;
typedef PDT::Locate_type                                    Locate_type;
typedef PDT::Iso_rectangle                                  Iso_rectangle;


int main()
{
  CGAL::Timer t;
  typedef CGAL::Creator_uniform_2<double, Point> Creator;
  CGAL::Random random(7);
  CGAL::Random_points_in_square_2<Point, Creator> in_square(.5, random);

  Iso_rectangle domain(0.0, 0.0, 1.0, 1.0); // The cube for the periodic domain
  int n = 4096;
  std::list<Point> pts;
  FILE *fp1;

  double x1, y1;
  PDT  PT2, PT3;

  fp1 = fopen("test.dat", "r");
  fprintf(stderr, " file opened %lf %lf \n", 1.0, 1.0);
  // Generating n random points
  for (int i = 0 ; i < n ; i++)
  {
      fscanf(fp1, "%lf %lf\n", &x1, &y1);
    /* Point p = *in_square; */
    /* in_square++; */
    pts.push_back(Point(x1, y1));
    /* fprintf(stderr, "%lf %lf \n", pts[i].x, pts[i].y); */
  }
  fclose(fp1);

  // Standard insertion
  t.start();
  /* for (int i = 0 ; i < n ; i++) */
    /* PT1.insert(pts[i]); */

  PDT PT1(pts.begin(), pts.end(), domain); // Put the domain with the constructor
  t.stop();
  std::cout << "  Time: " << t.time() << " sec. (Standard insertion)" << std::endl;
  t.reset();

  Face_iterator it = PT1.faces_begin(),
          beyond = PT1.faces_end();

  Vertex_iterator vt = PT1.vertices_begin(),
          vt_beyond = PT1.vertices_end();

  Face_handle face;
  Vertex_handle vtx;
  Locate_type lt;
  int li, i;

  fp1 = fopen("points__.dat", "w");
  while (vt != vt_beyond) {
      vtx = vt; // get face
    fprintf(fp1, "%lf %lf \n", PT1.point(vtx).x(), PT1.point(vtx).y());
     /* std::cout << a << "  "  << b << "  "  << c << std::endl; */
      ++vt; // advance the iterator
  }
  fclose(fp1);

  fp1 = fopen("simplices__.dat", "w");
  while (it != beyond) {
      face = it; // get face
      Point p, ptmp;
      int a=0, b=0, c = 0;
      double mod1, mod2, mod3, mod4;
      double xt1, yt1;
      vt = PT1.vertices_begin();
      i = 0;
      while (vt != vt_beyond) {
        vtx = vt; // get face
        p = PT1.point(vtx);
          mod1 = p.x()*p.x() + p.y()*p.y();

          xt1 = PT1.triangle(face)[0].x();
          yt1 = PT1.triangle(face)[0].y();
          if(xt1 > 1.0e0) xt1 = xt1 - 1.0e0;
          if(yt1 > 1.0e0) yt1 = yt1 - 1.0e0;
          mod2 = xt1*xt1 + yt1*yt1;

          xt1 = PT1.triangle(face)[1].x();
          yt1 = PT1.triangle(face)[1].y();
          if(xt1 > 1.0e0) xt1 = xt1 - 1.0e0;
          if(yt1 > 1.0e0) yt1 = yt1 - 1.0e0;
          mod3 = xt1*xt1 + yt1*yt1;

          xt1 = PT1.triangle(face)[2].x();
          yt1 = PT1.triangle(face)[2].y();
          if(xt1 > 1.0e0) xt1 = xt1 - 1.0e0;
          if(yt1 > 1.0e0) yt1 = yt1 - 1.0e0;
          mod4 = xt1*xt1 + yt1*yt1;

            if(fabs(mod1-mod2) < 1e-6) a = i;
            if(fabs(mod1-mod3) < 1e-6) b = i; 
            if(fabs(mod1-mod4) < 1e-6) c = i;      /* for (int i = 0 ; i < n ; i++) */
            ++vt;
            i = i + 1;
     } 
      fprintf(fp1, "%d %d %d \n", a, b, c);
      std::cout << a << "  "  << b << "  "  << c << std::endl;
      /* std::cout << PT1.triangle(face)  << std::endl; */
      ++it; // advance the iterator

  }

  fclose(fp1);

  // Iterator range insertion using spatial sorting but no dummy points
  /* t.start(); */
  /* PT2.insert(pts.begin(), pts.end()); // third parameter defaults to false */
  /* t.stop(); */
  /* std::cout << "  Time: " << t.time() << " sec. (with spatial sorting)" << std::endl; */
  /* t.reset(); */

  /* // Iterator range insertion using spatial sorting and dummy point heuristic */
  /* t.start(); */
  /* PT3.insert(pts.begin(), pts.end(), true); */
  /* t.stop(); */
  /* std::cout << "  Time: " << t.time() << " sec. (Dummy point heuristic)" << std::endl; */

  return 0;
}
