#include "curve.h"

/*!
\brief Creates a cubic curve.

\param x, y %Cubic parametric equations of the corresponding x, y coordinates.
*/
CubicCurve2::CubicCurve2(const Cubic& x, const Cubic& y)
{
    CubicCurve2::x = x;
    CubicCurve2::y = y;
}

/*!
\brief Computes the tangent to the cubic curve at a given position on the curve.

The tangent is the obtained by evaluating the first derivative.
\param t Position.
*/
Vector2 CubicCurve2::Tangent(const double& t) const
{
    return Vector2(x.Prime()(t), y.Prime()(t));
}

/*!
\brief Computes the normal vector to the cubic curve at a given position on the curve.

The normal is the obtained by evaluating the second derivative.
\param t Position.
*/
Vector2 CubicCurve2::Normal(const double& t) const
{
    return Vector2(x.Second()(t), y.Second()(t));
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector2 CubicCurve2::operator() (const double& t) const
{
    return Vector2(x(t), y(t));
}

/*!
\brief Evaluates curve point at a given location.
\param t Parameter.
*/
Vector2 CubicCurve2::Eval(const double& t) const
{
    return Vector2(x(t), y(t));
}

/*!
\brief Creates an Hermite cubic curve on interval [0,1] given
vertex locations and tangent vectors (in that order).
\param a,b End vertices of the Hermite cubic curve.
\param ta,tb Tangent vectors at the end vertices.
*/
CubicCurve2 CubicCurve2::Hermite(const Vector2& a, const Vector2& b, const Vector2& ta, const Vector2& tb)
{
    return CubicCurve2(Cubic::Hermite(a[0], b[0], ta[0], tb[0]), Cubic::Hermite(a[1], b[1], ta[1], tb[1]));
}

/*!
\brief Creates a Bezier cubic curve on interval [0,1] given
four vertex locations.
\param a,b,c,d Vertices of the Bezier curve.
*/
CubicCurve2 CubicCurve2::Bezier(const Vector2& a, const Vector2& b, const Vector2& c, const Vector2& d)
{
    return CubicCurve2(Cubic::Bezier(a[0], b[0], c[0], d[0]), Cubic::Bezier(a[1], b[1], c[1], d[1]));
}

/*!
\brief Compute the parameter corresponding to the argument curvilign absisca of a point on the curve.

\param s Argument curvilign absisca (should be between 0.0 and the length of the curve).
\param n Discretization of the integration interval.
*/
double CubicCurve2::U(const double& s, int n) const
{
    if (s == 0.0) return 0.0;
    const double epsilon = 1.0 / n;
    double ss = 0.0;
    Vector2 a = Vector2(x(0.0), y(0.0));
    for (double u = epsilon; u < 1.0; u += epsilon)
    {
        Vector2 b = Vector2(x(u), y(u));
        ss += Magnitude(b - a);
        if (ss > s)
        {
            return u - epsilon / 2.0;
        }
        a = b;
    }
    return 1.0;
}

/*!
\brief Compute the curvilign absisca of a point on the curve.

\param t Parameter that should be within [0,1] interval.
\param n Discretization of the integration interval.
*/
double CubicCurve2::S(const double& t, int n) const
{
    const double epsilon = 1.0 / n;
    double s = 0.0;
    Vector2 a = Vector2(x(0.0), y(0.0));
    for (double u = epsilon; u < t; u += epsilon)
    {
        Vector2 b = Vector2(x(u), y(u));
        s += Magnitude(b - a);
        a = b;
    }
    return s;
}
