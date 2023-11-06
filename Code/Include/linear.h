#pragma once

#include "basics.h"

class Linear {
protected:
  double c[2] = { 0.0 }; //!< Array of coefficients.
public:
  //! Empty.
  Linear() {}
  explicit Linear(const double&, const double&);
  explicit Linear(const double&);
  //! Empty.
  ~Linear() {}

  constexpr double& operator[] (int);
  constexpr double operator[] (int) const;

  // Unary operators
  Linear operator+ () const { return *this; }
  Linear operator- () const;

  Linear& operator+= (const Linear&);
  Linear& operator-= (const Linear&);
  Linear& operator*= (const double&);
  Linear& operator/= (const double&);

  // Binary operators
  friend Linear operator+ (const Linear&, const Linear&);
  friend Linear operator- (const Linear&, const Linear&);

  friend Linear operator* (const Linear&, const double&);
  friend Linear operator* (const double&, const Linear&);
  friend Linear operator/ (const Linear&, const double&);

  // Evaluates polynomial
  constexpr double operator()(const double&) const;
public:
  static const Linear Id; // Identity.
};

/*!
\brief Creates a linear function.

The function is defined as f(x) = a x + b.
\param a,b Coefficients.
*/
inline Linear::Linear(const double& a, const double& b)
{
  Linear::c[1] = a;
  Linear::c[0] = b;
}

/*!
\brief Creates a constant function.
\param a Constant.
*/
inline Linear::Linear(const double& a)
{
  Linear::c[1] = 0.0;
  Linear::c[0] = a;
}

//! Access class components.
inline constexpr double& Linear::operator[] (int i)
{
  return c[i];
}

//! Overloaded.
inline constexpr double Linear::operator[] (int i) const
{
  return c[i];
}

/*!
\brief Evaluates the linear function.
\param x Argument.
*/
inline constexpr double Linear::operator()(const double& x) const
{
  return c[0] + x * c[1];
}

//! Multiply a polynomial by a scalar value.
inline Linear operator* (const Linear& u, const double& e)
{
  return Linear(e * u[1], e * u[0]);
}

//! Overloaded.
inline Linear operator+ (const Linear& u, const Linear& v)
{
  return Linear(v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Linear operator- (const Linear& v, const Linear& u)
{
  return Linear(v[1] - u[1], v[0] - u[0]);
}

//! Destructive sum of two linear polynomials.
inline Linear& Linear::operator+=(const Linear& u)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] += u[i];
  }

  return *this;
}

//! Destructive difference of two linear polynomials.
inline Linear& Linear::operator-=(const Linear& u)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] -= u[i];
  }
  return *this;
}

//! Scale a linear polynomial by a double value.
inline Linear& Linear::operator*=(const double& e)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] *= e;
  }
  return *this;
}

//! Scales a linear polynomial by a double value.
inline Linear& Linear::operator/=(const double& e)
{
  for (int i = 0; i < 2; i++)
  {
    c[i] /= e;
  }
  return *this;
}

//! Unary.
inline Linear Linear::operator- () const
{
  return Linear(-c[1], -c[0]);
}

//! Overloaded.
inline Linear operator*(const double& a, const Linear& p)
{
  return p * a;
}

//! Overloaded.
inline Linear operator/(const Linear& p, const double& a)
{
  return p * (1.0 / a);
}
