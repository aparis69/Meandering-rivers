#pragma once

#include "basics.h"

class Quadric {
protected:
  double c[3] = { 0.0,0.0,0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Quadric() {}
  explicit Quadric(const double&, const double&, const double&);

  //! Empty.
  ~Quadric() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  // Unary operators
  Quadric operator+ () const { return *this; }
  Quadric operator- () const;

  Quadric& operator+= (const Quadric&);
  Quadric& operator-= (const Quadric&);
  Quadric& operator*= (const double&);
  Quadric& operator/= (const double&);

  // Binary operators
  friend Quadric operator+ (const Quadric&, const Quadric&);
  friend Quadric operator- (const Quadric&, const Quadric&);

  friend Quadric operator* (const Quadric&, const double&);
  friend Quadric operator* (const double&, const Quadric&);
  friend Quadric operator/ (const Quadric&, const double&);

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
  constexpr double Derivative(const double&) const;

  // Solve
  int Solve(double&, double&) const;
  int Solve(double*) const;
  int Solve(double*, const double&, const double&) const;

  void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;
  double Maximum(const double& = 0.0, const double& = 1.0) const;
  double Minimum(const double& = 0.0, const double& = 1.0) const;

  static Quadric Bezier(const double&, const double&, const double&);

public:
  static double epsilon; //!< Epsilon value used to check b<SUP>2</SUP>-4ac term in the root finding process.
};

//! Creates a quadric.
inline Quadric::Quadric(const double& a, const double& b, const double& c)
{
  Quadric::c[2] = a;
  Quadric::c[1] = b;
  Quadric::c[0] = c;
}

//! Access class components.
inline constexpr double& Quadric::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Quadric::operator[] (int i) const
{
  return c[i];
}

//! Evaluates the quadric. 
inline constexpr double Quadric::operator()(const double& x) const
{
  return c[0] + x * (c[1] + x * c[2]);
}

/*!
\brief Computes the derivative value at a given point.

If only one evaluation is to be made, this function is more efficient than Prime() and
then evaluating the derivative for a given input value.

\param x Real.

\sa Quadric::Prime
*/
inline constexpr double Quadric::Derivative(const double& x) const
{
  return 2.0 * c[2] * x + c[1];
}

/*!
\brief Destructive sum of two polynomials.
*/
inline Quadric& Quadric::operator+=(const Quadric& u)
{
  c[0] += u[0];
  c[1] += u[1];
  c[2] += u[2];

  return *this;
}

/*!
\brief Destructive difference of two polynomials.
*/
inline Quadric& Quadric::operator-=(const Quadric& u)
{
  c[0] -= u[0];
  c[1] -= u[1];
  c[2] -= u[2];

  return *this;
}

/*!
\brief Scale a polynomial by a double value.
*/
inline Quadric& Quadric::operator*=(const double& e)
{
  c[0] *= e;
  c[1] *= e;
  c[2] *= e;
  return *this;
}

/*!
\brief Scale a polynomial by a double value.
*/
inline Quadric& Quadric::operator/=(const double& e)
{
  c[0] /= e;
  c[1] /= e;
  c[2] /= e;
  return *this;
}

/*!
\brief Multiply a quadric by a scalar value.
*/
inline Quadric operator*(const Quadric& p, const double& e)
{
  return Quadric(e * p[2], e * p[1], e * p[0]);
}

//! Overloaded.
inline Quadric operator+(const Quadric& u, const Quadric& v)
{
  return Quadric(v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Quadric operator-(const Quadric& v, const Quadric& u)
{
  return Quadric(v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

//! Overloaded
inline Quadric Quadric::operator- () const
{
  return Quadric(-c[2], -c[1], -c[0]);
}

/*!
\brief Overloaded.
*/
inline Quadric operator*(const double& a, const Quadric& p)
{
  return Quadric(a * p[2], a * p[1], a * p[0]);
}

/*!
\brief Overloaded.
*/
inline Quadric operator/(const Quadric& p, const double& a)
{
  return p * (1.0 / a);
}
