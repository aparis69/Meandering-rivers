#pragma once

#include "basics.h"

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

    int CheckDegree() const;

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

    // Solve
    int Solve(double*) const;
    int Solve(double*, const double&, const double&) const;
    int SolveNormalized(double*) const;

    void Range(double&, double&, const double& = 0.0, const double& = 1.0) const;

    double K(const double&, const double&) const;

    static Cubic Hermite(const double&, const double&, const double&, const double&);
    static Cubic Bezier(const double&, const double&, const double&, const double&);
    static Cubic Spline(const double&, const double&, const double&, const double&);
    static double Interpolation(const double&, const double&, const double&, const double&, const double&);

    static double Smooth(const double&, const double&);
    static double SmoothCompact(const double&, const double&);
    static double Smooth(const double&);
    static double SmoothStep(const double&, const double&, const double&);

    static double Gaussian(const double&, const double&, const double&);
    static double GaussianThick(const double&, const double&, const double&, const double&);

    // Bernstein
    static Cubic Bernstein(int);
    static double Bernstein(int, const double&);
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
\brief Compute the value of a C<SUP>2</SUP> smooth interpolating function (1-x/r)<SUP>3</SUP>.

It is possible to implement Wyvill's smoothing kernel by passing the squared
distance to the skeleton and the squared radius to this function. For example,
given a box skeleton:
\code
Vector p(-2.0,3.0,5.0);
double d=Box(1.0).R(p); // Squared distance to a box
double r=4.0; // Radius;
double d=d>r*r?0.0:Cubic::Smooth(d,r*r); // Direct implementation of Wyvill's smoothing kernel
\endcode

\sa Smooth(const double&,const double&,const double&)

\param x Squared distance.
\param r Squared radius.
*/
inline double Cubic::Smooth(const double& x, const double& r)
{
    return (1.0 - x / r) * (1.0 - x / r) * (1.0 - x / r);
}

/*!
\brief Compactly supported smooth interpolating function.

\sa Cubic::Smooth(const double&, const double&), Quadric::SmoothCompact()

\param x Squared distance.
\param r Squared radius.
*/
inline double Cubic::SmoothCompact(const double& x, const double& r)
{
    return (x > r) ? 0.0 : (1.0 - x / r) * (1.0 - x / r) * (1.0 - x / r);
}

/*!
\brief Compute the value of a C<SUP>1</SUP> smooth interpolating function.

The cubic is defined as x<SUP>2</SUP>(3 - 2 x). Its first derivatives at 0.0 and 1.0 are 0.0.

The Lipschitz constant of the smooth quintic over [0,1] is &lambda;=3/2.

\param x Argument in [0,1].

\sa Quintic::Smooth(), Septic::Smooth()
*/
inline double Cubic::Smooth(const double& x)
{
    return x * x * (3.0 - 2.0 * x);
}

/*!
\brief Compute a smooth cubic step.

The code is slightly more efficient than:
\code
double y=Cubic::Smooth(Linear::Step(x,a,b));
\endcode
\param a, b Interval values.
\param x Input value.
\sa Quintic::SmoothStep(), Septic::SmoothStep()
*/
inline double Cubic::SmoothStep(const double& x, const double& a, const double& b)
{
    if (x < a)
    {
        return 0.0;
    }
    else if (x > b)
    {
        return 1.0;
    }
    else
    {
        return Cubic::Smooth((x - a) / (b - a));
    }
}

/*!
\brief Compute a compactly supported Gaussian-like pulse.

The function has C<SUP>1</SUP> continuity.
\sa Quintic::Gaussian()

\param c Center.
\param r Radius.
\param x Value.
*/
inline double Cubic::Gaussian(const double& c, const double& r, const double& x)
{
    double xc = fabs(x - c);
    if (xc > r) return 0.0;
    xc /= r;
    return Cubic::Smooth(1.0 - xc);
}

/*!
\brief Compute a compactly supported Gaussian-like pulse with a thick plateau.

The function has C<SUP>1</SUP> continuity.
\sa Quintic::GaussianThick()

The support has 2 ( r + t ) extent.

\param c Center.
\param t,r Thickness (plateau) and radius.
\param x Value.
*/
inline double Cubic::GaussianThick(const double& c, const double& t, const double& r, const double& x)
{
    double y = fabs(x - c);
    if (y > r + t) return 0.0;
    y -= t;
    if (y < 0.0) return 1.0;
    y /= r;
    return Cubic::Smooth(1.0 - y);
}
