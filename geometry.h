//
// Created by Nikita Mastinen on 15.11.2020.
//

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#endif //GEOMETRY_GEOMETRY_H

#include <utility>
#include <vector>
#include <cmath>

bool isDoublesEqual(const double a, const double b, double precision = 1e-11) {
  return -precision <= a - b && a - b <= precision;
}

struct Point;
struct Vector;

struct Point {
  double x = 0.0, y = 0.0;

  Point() = default;

  explicit Point(const Vector& a);

  Point(double x, double y): x(x), y(y) {};

  bool operator==(const Point& other) const {
    return isDoublesEqual(x, other.x) && isDoublesEqual(y, other.y);
  }

  bool operator!=(const Point& other) const {
    return !isDoublesEqual(x, other.x) || !isDoublesEqual(y, other.y);
  }
};

struct Vector {
  double x = 0.0, y = 0.0;

  Vector() = default;

  Vector(double x, double y): x(x), y(y) {};

  Vector(Point a, Point b): x(b.x - a.x), y(b.y - a.y) {};

  Vector& normalize() {
    double k = std::pow(x * x + y * y, 0.5);
    x /= k, y /= k;
    return *this;
  }

  bool operator==(const Vector& other) const {
    return isDoublesEqual(x, other.x) && isDoublesEqual(y, other.y);
  }

  bool operator!=(const Vector& other) const {
    return !isDoublesEqual(x, other.x) || !isDoublesEqual(y, other.y);
  }
};

Vector operator*(const Vector& a, const double k) {
  return {a.x * k, a.y * k};
}

Vector operator+(const Vector& a, const Vector& b) {
  return {a.x + b.x, a.y + b.y};
}

double operator*(const Vector& a, const Vector& b) {
  return a.x * b.y - b.x * a.y;
}

Point::Point(const Vector& a): x(a.x), y(a.y) {}

class Line {
public:
  double a = 1.0, b = 0.0, c = 0.0;

  Line() = default;

  Line(const double a, const double b, const double c): a(a), b(b), c(c) {};

  Line(const Point& u, const Point& v): a(v.y - u.y),
      b(u.x - v.x),
      c(-u.x * v.y + v.x * u.y) {};

  Line(const double k, const double b): a(k), b(-1.0), c(b) {};

  Line(const Point& a, const double & k): Line(k, a.y - a.x * k) {};
};

Point intersect(Line l, Line m) {
  double det = l.a * m.b - m.a * l.b;
  return {(m.c * l.b - l.c * m.b) / det, (-m.c * l.a + l.c * m.a) / det};
}

class Shape {
  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool operator==(const Shape& another) const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;

  virtual Shape& rotate(const Point& center, double angle) = 0;

  virtual Shape& reflex(const Point& center) = 0;

  virtual Shape& reflex(Line axis) = 0;

  virtual Shape& scale(Point center, double coefficient) = 0;
};

class Polygon: public Shape {
protected:
  std::vector<Point> vertices;
public:
  Polygon() = default;

  Polygon(const std::vector<Point>& vertices): vertices(vertices) {};

  Polygon(const std::initializer_list<Point>& l) : vertices(l) {};

  size_t verticesCount() const {
    return vertices.size();
  }

  const std::vector<Point>& getVertices() const {
    return vertices;
  }

  bool isConvex() const {
   //TODO implement!!!
  }

  double perimeter() const override {
      //TODO implement!!!
  }

  double area() const override {
      //TODO implement!!!
  }

  bool operator==(const Shape& another) const override {
      //TODO implement!!!
  }

  bool isCongruentTo(const Shape& another) const override {
      //TODO implement!!!
  }

  bool isSimilarTo(const Shape& another) const override {
    //TODO implement!!!
  }

  bool containsPoint(const Point& point) const override {
    //TODO implement!!!
  }

  Shape& rotate(const Point& center, double angle) override {
    //TODO implement!!!
  }

  Shape& reflex(const Point& center) override {
    //TODO implement!!!
  }

  Shape& reflex(Line axis) override {
    //TODO implement!!!
  }

  Shape& scale(Point center, double coefficient) override {
    //TODO implement!!!
  }
};

class Ellipse: public Shape {
protected:
  std::pair<Point, Point> focus;
  double sum_of_distances = 0.0;
public:
  Ellipse() = default;

  Ellipse(const std::pair<Point, Point>& focus, const double sum_of_distances):
    focus(focus), sum_of_distances(sum_of_distances) {};

  const std::pair<Point,Point>& focuses() const {
    return focus;
  }

  std::pair<Line, Line> directrices() const {
    //TODO implement!!!
  }

  Point eccentricity() const {
    //TODO implement!!!
  }

  Point center() const {
    //TODO implement!!!
  }

  bool isConvex() const {
    //TODO implement!!!
  }

  double perimeter() const override {
    //TODO implement!!!
  }

  double area() const override {
    //TODO implement!!!
  }

  bool operator==(const Shape& another) const override {
    //TODO implement!!!
  }

  bool isCongruentTo(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool isSimilarTo(const Shape& another) const override {
      //TODO implement!!!
  }

  bool containsPoint(const Point& point) const override {
      //TODO implement!!!
  }

  Shape& rotate(const Point& center, double angle) override {
    //TODO implement!!!
  }

  Shape& reflex(const Point& center) override {
    //TODO implement!!!
  }

  Shape& reflex(Line axis) override {
    //TODO implement!!!
  }

  Shape& scale(Point center, double coefficient) override {
    //TODO implement!!!
  }
};

class Circle: public Ellipse {
  Point center;
  double r = 0.0;
public:
  Circle() = default;

  Circle(Point center, double r): Ellipse({center, center}, 2.0 * r), center(center), r(r) {};

  double radius() const {
    return r;
  }

  double perimeter() const override {
    //TODO implement!!!
  }

  double area() const override {
    //TODO implement!!!
  }

  bool operator==(const Shape& another) const override {
    //TODO implement!!!
  }

  bool isCongruentTo(const Shape& another) const override {
    //TODO implement!!!
  }

  bool isSimilarTo(const Shape& another) const override {
    //TODO implement!!!
  }

  bool containsPoint(const Point& point) const override {
    //TODO implement!!!
  }
};

class Rectangle: public Polygon {
public:
  Rectangle(Point a, Point c, double k) {
    //TODO Implement!!!
  };

  Point center() {
    return Point((vertices[0].x + vertices[2].x) / 2, (vertices[0].y + vertices[2].y) / 2);
  }

  std::pair<Line, Line> diagonals() {
    return {Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3])};
  }
};

class Square: public Rectangle {
public:
  Square(Point a, Point c): Rectangle(a, c, 1.0) {};

  Circle circumscribedCircle() {
    return Circle(center(), std::abs(vertices[0].x - vertices[1].x) / std::pow(2.0, 0.5));
  }

  Circle inscribedCircle() {
    return Circle(center(), std::abs(vertices[0].x - vertices[1].x) / 2);
  }
};

class Triangle: public Polygon {
public:
  Circle circumscribedCircle() {
    //TODO
  }

  Circle inscribedCircle() {
    //TODO
  }

  Point centroid() {
    //TODO
  }

  Point orthocenter() {
    //TODO
  }

  Line EulerLine() {
    //TODO
  }

  Circle ninePointsCircle() {
    //TODO
  }
};




