#pragma once

#include "basics.h"
#include "quadric.h"
#include "linear.h"

class Cubic{
protected:
  double c[4] = { 0.0,0.0,0.0,0.0 }; //!< Array of coefficients.
public:
    //! Empty.
    Cubic() {}
    explicit Cubic(const double&, const double&, const double&, const double&);
    //! Empty.
    ~Cubic() {}

    constexpr double& operator[] (int);
    constexpr double operator[] (int) const;

    // Unary operators

    // Overloaded out of consistency.
    Cubic operator+ () const { return *this; }
    Cubic operator- () const;

    Cubic& operator+= (const Cubic&);
    Cubic& operator-= (const Cubic&);
    Cubic& operator*= (const double&);
    Cubic& operator/= (const double&);

    // Binary operators
    friend Cubic operator+ (const Cubic&, const Cubic&);
    friend Cubic operator- (const Cubic&, const Cubic&);

    friend Cubic operator* (const Cubic&, const double&);
    friend Cubic operator* (const double&, const Cubic&);
    friend Cubic operator/ (const Cubic&, const double&);

    // Evaluates cubic
    constexpr double operator()(const double&) const;
    constexpr double Derivative(const double&) const;

    Quadric Prime() const;
    Linear Second() const;

    // Solve
    int Solve(double*) const;
    int Solve(double*, const double&, const double&) const;
    int SolveNormalized(double*) const;

    static Cubic Hermite(const double&, const double&, const double&, const double&);
    static Cubic Bezier(const double&, const double&, const double&, const double&);
  private:
    static double epsilon; //!< \htmlonly\epsilon;\endhtmlonly value used to check if the cubic term is so close to 0 that the cubic is in fact a quadric.
};

/*!
\brief Creates a cubic.

Coefficients start from the highest degree to the lowest.
\param a, b, c, d Coefficients of the cubic a x<SUP>3</SUP>+b x<SUP>2</SUP>+c x+d.
*/
inline Cubic::Cubic(const double& a, const double& b, const double& c, const double& d)
{
    Cubic::c[3] = a;
    Cubic::c[2] = b;
    Cubic::c[1] = c;
    Cubic::c[0] = d;
}

//! Access class components.
inline constexpr double& Cubic::operator[] (int i)
{
    return c[i];
}

//! Overloaded.
inline constexpr double Cubic::operator[] (int i) const
{
    return c[i];
}

/*!
\brief Evaluates the cubic.

\param x Value.
*/
inline constexpr double Cubic::operator()(const double& x) const
{
    return c[0] + x * (c[1] + x * (c[2] + x * c[3]));
}

/*!
\brief Computes the derivative of the cubic.

This function is more efficient than using Prime() and
then evaluating the derivative for a given input value.
\param x Real.
*/
inline constexpr double Cubic::Derivative(const double& x) const
{
    return c[1] + x * (2.0 * c[2] + x * 3.0 * c[3]);
}

//! Overloaded.
inline Cubic operator+ (const Cubic& u, const Cubic& v)
{
    return Cubic(v[3] + u[3], v[2] + u[2], v[1] + u[1], v[0] + u[0]);
}

//! Overloaded.
inline Cubic operator- (const Cubic& v, const Cubic& u)
{
    return Cubic(v[3] - u[3], v[2] - u[2], v[1] - u[1], v[0] - u[0]);
}

//! Unary.
inline Cubic Cubic::operator- () const
{
    return Cubic(-c[3], -c[2], -c[1], -c[0]);
}

//! Multiply a cubic by a scalar value.
inline Cubic operator* (const Cubic& u, const double& e)
{
    return Cubic(e * u[3], e * u[2], e * u[1], e * u[0]);
}

//! Overloaded.
inline Cubic operator*(const double& a, const Cubic& p)
{
    return p * a;
}

//! Overloaded.
inline Cubic operator/(const Cubic& p, const double& a)
{
    return p * (1.0 / a);
}

/*!
\brief Destructive sum of two cubics.
*/
inline Cubic& Cubic::operator+=(const Cubic& u)
{
    for (int i = 0; i < 4; i++)
    {
        c[i] += u[i];
    }

    return *this;
}

/*!
\brief Destructive difference of two cubics.
*/
inline Cubic& Cubic::operator-=(const Cubic& u)
{
    for (int i = 0; i < 4; i++)
    {
        c[i] -= u[i];
    }
    return *this;
}

/*!
\brief Scale a cubic by a double value.
*/
inline Cubic& Cubic::operator*=(const double& e)
{
    for (int i = 0; i < 4; i++)
    {
        c[i] *= e;
    }
    return *this;
}

/*!
\brief Scale a cubic by a double value.
*/
inline Cubic& Cubic::operator/=(const double& e)
{
    for (int i = 0; i < 4; i++)
    {
        c[i] /= e;
    }
    return *this;
}

/*!
\brief Computes the first derivative of a cubic, which is a quadric.
*/
inline Quadric Cubic::Prime() const
{
  return Quadric(3.0 * c[3], 2.0 * c[2], c[1]);
}

/*!
\brief Computes the second derivative of a cubic, which is a linear polynomial.
*/
inline Linear Cubic::Second() const
{
  return Linear(6.0 * c[3], 2.0 * c[2]);
}
