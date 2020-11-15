//
// Created by Nikita Mastinen on 15.11.2020.
//

#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

#endif //GEOMETRY_GEOMETRY_H

#include <utility>
#include <vector>

bool isDoublesEqual(const double a, const double b, double precision = 1e-11) {
  return -precision <= a - b && a - b <= precision;
}

struct Point {
  double x = 0.0, y = 0.0;

  Point() = default;

  Point(int x, int y): x(x), y(y) {};

  bool operator==(const Point& other) const {
    return isDoublesEqual(x, other.x) && isDoublesEqual(y, other.y);
  }

  bool operator!=(const Point& other) const {
    return !isDoublesEqual(x, other.x) || !isDoublesEqual(y, other.y);
  }
};

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

class Shape {
  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool operator==(const Shape& another) const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;
};

class Polygon: public Shape {
private:
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

  virtual double perimeter() const override {
      //TODO implement!!!
  }

  virtual double area() const override {
      //TODO implement!!!
  }

  virtual bool operator==(const Shape& another) const override {
      //TODO implement!!!
  }

  virtual bool isCongruentTo(const Shape& another) const override {
      //TODO implement!!!
  }

  virtual bool isSimilarTo(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool containsPoint(const Point& point) const override {
    //TODO implement!!!
  }
};

class Ellipse: public Shape {
  std::pair<Point, Point> focus;
  double sum_of_distances = 0;
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

  virtual double perimeter() const override {
    //TODO implement!!!
  }

  virtual double area() const override {
    //TODO implement!!!
  }

  virtual bool operator==(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool isCongruentTo(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool isSimilarTo(const Shape& another) const override {
      //TODO implement!!!
  }

  virtual bool containsPoint(const Point& point) const override {
      //TODO implement!!!
  }
};

class Circle: public Ellipse {
//TODO implement!
public:

  virtual double perimeter() const override {
    //TODO implement!!!
  }

  virtual double area() const override {
    //TODO implement!!!
  }

  virtual bool operator==(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool isCongruentTo(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool isSimilarTo(const Shape& another) const override {
    //TODO implement!!!
  }

  virtual bool containsPoint(const Point& point) const override {
    //TODO implement!!!
  }
};




