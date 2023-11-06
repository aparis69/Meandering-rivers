#include "cubic.h"

double Cubic::epsilon = 1.0e-10;

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
