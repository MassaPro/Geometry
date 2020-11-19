#include <iostream>
#include <vector>
#include "geometry.h"

using namespace std;

int main() {
  Polygon ppp = Polygon({Point(18.1565, -1.82785), Point(4.46915, 5.00267), Point(-4.21555, 7.36405), Point(-4.46915, -5.00267), Point(2.10777, -3.68202), Point(7.89757, -5.25627)});
  Polygon ddd = Polygon({Point(-2, 2),  Point(1, 2),  Point(6, 1),  Point(3, -1), Point(1,  -1), Point( -1, -2)});

  Polygon a = Rectangle(Point(1, 1), Point(2,2), 1);
  Polygon b = Rectangle(Point(0, 0), Point(0, 2), 1);

  //cout << "is_similar " << a.containsPoint({1.5, 1.5});
  //cout << (Point(1, 2) == Point(1, 2));
  Triangle t = Triangle({1, 1}, {5, 5}, {5, 3});
  print(t.inscribedCircle().center());
  return 0;
}
