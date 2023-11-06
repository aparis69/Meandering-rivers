#pragma once

#include "basics.h"
#include "cubic.h"

class CubicCurve2
{
protected:
    Cubic x, y; //!< %Cubic polynomial functions for every coordinate.
public:
    //! Empty.
    CubicCurve2() {}
    explicit CubicCurve2(const Cubic&, const Cubic&);

    //! Empty.
    ~CubicCurve2() {}

    // Access class components
    Cubic& operator[] (int);
    Cubic operator[] (int) const;

    // Computes the tangent curve equation
    Vector2 Tangent(const double&) const;
    Vector2 Normal(const double&) const;

    // Compute the curvilign abscissa of a point on the curve
    double S(const double&, int = 256) const;
    double U(const double&, int = 256) const;

    // Evaluates curve
    Vector2 operator()(const double&) const;
    Vector2 Eval(const double&) const;

public:
    static CubicCurve2 Hermite(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
    static CubicCurve2 Bezier(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
};

inline Cubic& CubicCurve2::operator[] (int i)
{
    if (i == 0) { return x; }
    else { return y; }
}

inline Cubic CubicCurve2::operator[] (int i) const
{
    if (i == 0) { return x; }
    else { return y; }
}



class CubicCurve2Set
{
protected:
  std::vector<CubicCurve2> curve; //!< Set of 2D cubic curves.
  std::vector<double> lengths; //!< Length of every curve.
  double length = 0.0; //!< Total length.
public:
  CubicCurve2Set();
  explicit CubicCurve2Set(const std::vector<CubicCurve2>&);
  explicit CubicCurve2Set(const std::vector<Vector2>&, const Vector2& = Vector2(0), const Vector2 & = Vector2(0));

  //! Empty
  ~CubicCurve2Set() {}

  // Access to members
  int Size() const;
  CubicCurve2 operator()(int) const;
  std::vector<Vector2> GetDiscretisation(const double&) const;
  std::vector<Vector2> GetDiscretisation(const double&, std::vector<Vector2>&) const;
  int U(const double&, double&) const;
};

inline CubicCurve2 CubicCurve2Set::operator()(int i) const
{
  return curve.at(i);
}

inline int CubicCurve2Set::Size() const
{
  return int(curve.size());
}
