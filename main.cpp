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

  cout << c.getVertices()[0].x;
  return 0;
}
