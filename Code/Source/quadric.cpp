#include "quadric.h"

double Quadric::epsilon = 1.0e-16;

/*!
\brief This function computes the range of values taken by a quadric over a given interval.

Bascially computes the roots of the first derivative, and evaluates the polynomial
at the roots if they are within the interval bounds.

\param a,b Interval.
\param x,y Returned range.
*/
void Quadric::Range(double& x, double& y, const double& a, const double& b) const
{
  x = (*this)(a);
  y = (*this)(b);

  if (x > y)
    std::swap(x, y);

  double t = -c[1] / (c[2] + c[2]);

  if ((t > a) && (t < b))
  {
    double s = (*this)(t);
    if (s < x)
    {
      x = s;
    }
    else if (s > y)
    {
      y = s;
    }
  }
}

/*!
\brief This function computes the maximum value taken by a quadric over a given interval.

\sa Range
\param a,b Interval.
*/
double Quadric::Maximum(const double& a, const double& b) const
{
  double x = Math::Max((*this)(a), (*this)(b));

  double t = -c[1] / (c[2] + c[2]);

  if ((t > a) && (t < b))
  {
    double s = (*this)(t);
    if (s > x)
    {
      x = s;
    }
  }
  return x;
}

/*!
\brief This function computes the minimum value taken by a quadric over a given interval.

\sa Range
\param a,b Interval.
*/
double Quadric::Minimum(const double& a, const double& b) const
{
  double x = Math::Min((*this)(a), (*this)(b));

  double t = -c[1] / (c[2] + c[2]);

  if ((t > a) && (t < b))
  {
    double s = (*this)(t);
    if (s < x)
    {
      x = s;
    }
  }
  return x;
}

/*!
\brief Solve the quadratic equation over a given interval.

This function store the
sorted roots in an array and returns the number of roots.

\param roots The array of roots.
\param a,b The interval.
*/
int Quadric::Solve(double* roots, const double& a, const double& b) const
{
  int r = Solve(roots);
  int j = 0;

  for (int i = 0; i < r; i++)
  {
    if ((roots[i] > a) && (roots[i] < b))
    {
      roots[j] = roots[i];
      j++;
    }
  }

  return j;
}

/*!
\brief Solve quadratic equations. This function store the
sorted roots in an array and returns the number of roots.

\param y The array of roots.
*/
int Quadric::Solve(double* y) const
{
  double a = Quadric::c[2];
  double b = -Quadric::c[1];
  double c = Quadric::c[0];
  if (a == 0.0)
  {
    if (b == 0.0)
      return 0;
    y[0] = y[1] = c / b;
    return 1;
  }
  double d = b * b - 4.0 * a * c;
  // No roots
  if (d < 0.0)
  {
    return 0;
  }
  // One root
  else if (fabs(d) < epsilon)
  {
    y[0] = y[1] = 0.5 * b / a;
    return 1;
  }

  d = sqrt(d);
  double t = 0.5 / a;
  if (t > 0.0)
  {
    y[0] = (b - d) * t;
    y[1] = (b + d) * t;
  }
  else
  {
    y[0] = (b + d) * t;
    y[1] = (b - d) * t;
  }
  return 2;
}

/*!
\brief Solve quadratic equations.

This function uses two
double arguments to store the roots and returns the number of roots.

\param u, v The two of roots, if they exist.
*/
int Quadric::Solve(double& u, double& v) const
{
  double a = Quadric::c[2];
  double b = -Quadric::c[1];
  double c = Quadric::c[0];

  if (a == 0.0)
  {
    if (b == 0.0)
      return 0;
    u = v = c / b;
    return 1;
  }
  double d = b * b - 4.0 * a * c;
  if (d < 0.0)
    return 0;
  else if (fabs(d) < epsilon)
  {
    u = v = 0.5 * b / a;
    return 1;
  }
  d = sqrt(d);
  double t = 0.5 / a;
  if (t > 0.0)
  {
    u = (b - d) * t;
    v = (b + d) * t;
  }
  else
  {
    u = (b + d) * t;
    v = (b - d) * t;
  }
  return 2;
}

/*!
\brief Creates a quadratic Bezier polynomial on interval [0,1].
\param a,b,c The parameters of the Bézier polynomial.

Parameters are as follows:

c(t)=a (1-t)<SUP>2</SUP>+2b t(1-t)+c t<SUP>2</SUP>=(a-2b+c) t<SUP>2</SUP>+2(b-a) t+a
*/
Quadric Quadric::Bezier(const double& a, const double& b, const double& c)
{
  return Quadric(a - 2.0 * b + c, 2.0 * (b - a), a);
}
