#include "cubic.h"

double Cubic::epsilon = 1.0e-10;

/*!
\brief Check the degree of the cubic.
*/
int Cubic::CheckDegree() const
{
    int n = 3;
    while (n > 0)
    {
        if (c[n] != 0.0)
            return n;
        n--;
    }
    return n;
}

/*!
\brief Search the roots of a polynomial equation over a given interval.

\param roots Array for storing the roots.
\param a,b Interval limits.
*/
int Cubic::Solve(double* roots, const double& a, const double& b) const
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
\brief Solve the cubic equation.

This function stores the sorted roots in an array and returns the number of roots.

\param y Array for storing the roots.
*/
int Cubic::Solve(double* y) const
{
    double a1, a2, a3;

    double a0 = c[3];

    if (fabs(a0) < epsilon)
    {
        return Quadric(c[2], c[1], c[0]).Solve(y);
    }
    else
    {
        if (a0 != 1.0)
        {
            a1 = c[2] / a0;
            a2 = c[1] / a0;
            a3 = c[0] / a0;
        }
        else
        {
            a1 = c[2];
            a2 = c[1];
            a3 = c[0];
        }
    }

    double A2 = a1 * a1;
    double Q = (A2 - 3.0 * a2) / 9.0;
    double R = (a1 * (A2 - 4.5 * a2) + 13.5 * a3) / 27.0;
    double Q3 = Q * Q * Q;
    double R2 = R * R;
    double d = Q3 - R2;
    double an = a1 / 3.0;

    if (d >= 0.0)
    {
        // Three double roots

        d = R / sqrt(Q3);

        double theta = acos(d) / 3.0;
        double sQ = -2.0 * sqrt(Q);

        y[0] = sQ * cos(theta) - an;
        y[1] = sQ * cos(theta + 2.0943951023931954923084) - an; // 2.0*Math::Pi/3.0
        y[2] = sQ * cos(theta + 4.1887902047863909846168) - an; // 4.0*Math::Pi/3.0

        return 3;
    }
    else
    {
        double sQ = pow(sqrt(R2 - Q3) + fabs(R), 1.0 / 3.0);

        if (R < 0)
        {
            y[0] = (sQ + Q / sQ) - an;
        }
        else
        {
            y[0] = -(sQ + Q / sQ) - an;
        }

        return 1;
    }
}

/*!
\brief Overloaded output-stream operator.
\param s Stream.
\param c %Cubic.
*/
std::ostream& operator<<(std::ostream& s, const Cubic& c)
{
    s << "Cubic(" << c[3] << ',' << c[2] << ',' << c[1] << ',' << c[0] << ')';
    return s;
}

/*!
\brief Creates an Hermite cubic polynomial.

Interval parameterization is unit.

\param a, b Values for t=0 and t=1.
\param ta, tb Derivative values for t=0 and t=1.

\sa Math::Cubic()
*/
Cubic Cubic::Hermite(const double& a, const double& b, const double& ta, const double& tb)
{
    return Cubic(tb + ta + 2.0 * (a - b), -tb - 2.0 * ta + 3.0 * (b - a), ta, a);
}

/*!
\brief Creates a cubic Bezier polynomial.

Interval parameterization is unit.

\param a, b, c, d Control values.
*/
Cubic Cubic::Bezier(const double& a, const double& b, const double& c, const double& d)
{
    return Cubic(d - a + 3.0 * (b - c), 3.0 * (a + c) - 6.0 * b, 3.0 * (b - a), a);
}

/*!
\brief Creates a cubic Spline polynomial.

Interval parameterization is unit.

\param a, b, c, d Control values.
*/
Cubic Cubic::Spline(const double& a, const double& b, const double& c, const double& d)
{
    return Cubic((-a + 3.0 * b - 3.0 * c + d) / 6.0, (3.0 * a - 6.0 * b + 3.0 * c) / 6.0, (-3.0 * a + 3.0 * c) / 6.0, (a + 4.0 * b + c) / 6.0);
}

/*!
\brief Bicubic interpolation based on four values.
\param a,b,c,d Four input values.
\param x Interpolant.
\author Lois Paulin
*/
double Cubic::Interpolation(const double& x, const double& a, const double& b, const double& c, const double& d)
{
    double a_m1 = ((-0.5 * x + 1.0) * x - 0.5) * x;
    double a_0 = (1.5 * x - 2.5) * x * x + 1.0;
    double a_1 = ((-1.5 * x + 2.0) * x + 0.5) * x;
    double a_2 = (0.5 * x - 0.5) * x * x;

    return a_m1 * a + a_0 * b + a_1 * c + a_2 * d;
}

/*!
\brief Compute the Lipschitz constant of the cubic.

\param a,b Interval.
*/
double Cubic::K(const double& a, const double& b) const
{
    double x, y;

    Quadric quadric = Prime();
    quadric.Range(x, y, a, b);
    return Math::Max(fabs(x), fabs(y));
}

/*!
\brief Compute the range of values taken by the cubic over a given interval.

Computes the roots of the first derivative, and evaluates the polynomial
at the roots if they are within the interval bounds.
\param a,b Interval.
\param x,y Returned range.
*/
void Cubic::Range(double& x, double& y, const double& a, const double& b) const
{
    double r[2];

    // Self refence
    const Cubic& cubic = *this;
    x = cubic(a);
    y = cubic(b);

    Math::Sort(x, y);

    // Compute derivative
    Quadric p = Prime();

    // Find roots
    int n = p.Solve(r);

    for (int i = 0; i < n; i++)
    {
        if ((r[i] > a) && (r[i] < b))
        {
            double s = cubic(r[i]);
            if (s < x)
            {
                x = s;
            }
            if (s > y)
            {
                y = s;
            }
        }
    }
}

/*!
\brief Solve the normalized cubic equation, i.e., the highest coefficient is 1.0.

This function stores the sorted roots in an array and returns the number of roots.

\param y The array of roots.
*/
int Cubic::SolveNormalized(double* y) const
{
    double A2 = c[2] * c[2];
    double Q = (A2 - 3.0 * c[1]) / 9.0;
    double R = (c[2] * (A2 - 4.5 * c[1]) + 13.5 * c[0]) / 27.0;
    double Q3 = Q * Q * Q;
    double R2 = R * R;
    double d = Q3 - R2;
    double an = c[2] / 3.0;

    if (d >= 0.0)
    {
        // Three double roots
        d = R / sqrt(Q3);

        double theta = acos(d) / 3.0;
        double sQ = -2.0 * sqrt(Q);

        y[0] = sQ * cos(theta) - an;
        y[1] = sQ * cos(theta + Math::TwoPiOverThree) - an;
        y[2] = sQ * cos(theta + Math::FourPiOverThree) - an;

        return 3;
    }
    else
    {
        double sQ = pow(sqrt(R2 - Q3) + fabs(R), 1.0 / 3.0);

        if (R < 0)
        {
            y[0] = (sQ + Q / sQ) - an;
        }
        else
        {
            y[0] = -(sQ + Q / sQ) - an;
        }

        return 1;
    }
}

/*!
\brief Creates a quadric Bernstein polynomial.
\param p p-th polynomial.
*/
Cubic Cubic::Bernstein(int p)
{
    if (p == 0)
        return Cubic(1.0, -3.0, 3.0, -1.0);
    else if (p == 1)
        return Cubic(3.0, -6.0, 3.0, 0.0);
    else if (p == 2)
        return Cubic(-3.0, 3.0, 0.0, 0.0);
    else // if (p == 3)
        return Cubic(-3.0, 3.0, 0.0, 0.0);
}

/*!
\brief Compute the quadric Bernstein polynomial for a given value.
\param p p-th polynomial.
\param x Real.
*/
double Cubic::Bernstein(int p, const double& x)
{
    if (p == 3)
    {
        return x * x * x;
    }
    else
    {
        double t = 1.0 - x;
        if (p == 0)
        {
            return t * t * t;
        }
        else if (p == 1)
        {
            return 3.0 * x * t * t;
        }
        else
        {
            return 3.0 * x * x * t;
        }
    }
}
