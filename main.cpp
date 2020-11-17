#include <iostream>
#include <vector>
#include "geometry.h"

using namespace std;
int main() {
  Line l(Point(1, 2), 6);
  std::cout << l.a << ' ' << l.b << ' ' << l.c;
  Point a;
  Polygon c;
  Circle u;
  Ellipse h = u;

  Vector uu = Vector(Point(1, 1), Point());

  auto cc = Point(uu);
  uu.normalize();
  cout << (uu * 2).x << ' ' << (uu * 2).y;
  return 0;
}
