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

    // Computes the length of the cubic curve
    double Length(int = 256) const;
    double Length(const double&, const double&, int = 256) const;

    // Evaluates curve
    Vector2 operator()(const double&) const;
    Vector2 Eval(const double&) const;

    // Curvature 
    double Curvature(const double&) const;
    double Curvature(const double&, const double&) const;

    // Get bounding box
    Box2D GetBox() const;

    // Compute distance between point and curve
    double R(const Vector2&, double&) const;

    Vector2 BezierControl(int) const;

public:
    static CubicCurve2 Hermite(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
    static CubicCurve2 Bezier(const Vector2&, const Vector2&, const Vector2&, const Vector2&);
};

//! Access curve components
inline Cubic& CubicCurve2::operator[] (int i)
{
    if (i == 0) { return x; }
    else { return y; }
}

//! Access curve components
inline Cubic CubicCurve2::operator[] (int i) const
{
    if (i == 0) { return x; }
    else { return y; }
}

/*!
\brief Compute the i-th Bézier control point of the curve.

Computations have been optimized.
\param i Index.
*/
inline Vector2 CubicCurve2::BezierControl(int i) const
{
    if (i == 0)
    {
        // Equivalent to c(0.0)
        return Vector2(x[0], y[0]);
    }
    else if (i == 1)
    {
        return Vector2(x[0] + x[1] / 3.0, y[0] + y[1] / 3.0);
    }
    else if (i == 2)
    {
        return Vector2(x[0] + (2.0 * x[1] + x[2]) / 3.0, y[0] + (2.0 * y[1] + y[2]) / 3.0);
    }
    else
    {
        // Equivalent to c(1.0)
        return Vector2(x[0] + x[1] + x[2] + x[3], y[0] + y[1] + y[2] + y[3]);
    }
}