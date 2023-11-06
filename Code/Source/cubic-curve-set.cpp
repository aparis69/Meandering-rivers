#include "curve.h"

/*!
\class CubicCurve2Set curveset.h
\brief Piecewise cubic curves in the plane.

The length of every curve CubicCurve2Set::lengths and the total length CubicCurve2Set::length are computed
to speed up computations (internal optimization).
\ingroup PlanarGroup
*/

/*!
\brief Creates an empty piecewise quadric curve.
*/
CubicCurve2Set::CubicCurve2Set()
{
}

/*!
\brief Creates the piecewise quadric curve given a set of CubicCurve2.
\param control Set of cubic curves.
*/
CubicCurve2Set::CubicCurve2Set(const std::vector<CubicCurve2>& control)
{
  curve.clear();
  curve = control;

  lengths.clear();
  lengths.resize(curve.size());
  length = 0.0;
  for (int i = 0; i < curve.size(); i++)
  {
    lengths[i] = curve[i].S(1.0, 100);
    length += lengths[i];
  }
}

/*!
\brief Creates the piecewise cubic curve given a set of control points with the Cattmul-Rom construction.
\param control Set of 2D points.
\param ta Starting tangent vector.
\param tb Ending tangent vector.
*/
CubicCurve2Set::CubicCurve2Set(const std::vector<Vector2>& control, const Vector2& ta, const Vector2& tb)
{
  if (control.size() >= 3)
  {
    for (int i = 0; i < control.size() - 1; i++)
    {
      Vector2 t1;
      Vector2 t2;
      if (i == 0)
      {
        t1 = ta;
      }
      else
      {
        t1 = 0.5 * (control[i + 1] - control[i - 1]);
      }
      if (i == control.size() - 2)
      {
        t2 = tb;
      }
      else
      {
        t2 = 0.5 * (control[i + 2] - control[i]);
      }
      curve.push_back(CubicCurve2::Hermite(control[i], control[i + 1], t1, t2));
    }
    lengths.clear();
    lengths.resize(curve.size());
    length = 0.0;
    for (int i = 0; i < curve.size(); i++)
    {
      lengths[i] = curve[i].S(1.0, 100);
      length += lengths[i];
    }
  }
}

/*!
\brief Generates a discretization of the curve with a linear curvilign absisca parameterization.

\param step Stepping distance.
\param tangents Tangents to the curve at sample points.
*/
std::vector<Vector2> CubicCurve2Set::GetDiscretisation(const double& step, std::vector<Vector2>& tangents) const
{
  double niveau = 0.0;
  int idCurve = 0;

  // Set of points
  std::vector<Vector2> p;
  Vector2 last;

  // Start with the very first vertex of the piecewise curve
  p.push_back(curve[0](0));
  tangents.push_back(curve[0].Tangent(0.0001));
  last = curve[0](0);

  while (true)
  {
    niveau += step; // Position suivante sur la trajectoire

    double u;
    idCurve = U(niveau, u);

    if (idCurve != -1)
    {
      if (last != curve[idCurve](u))
      {
        p.push_back(curve[idCurve](u));  // Condition d'arret
        tangents.push_back(curve[idCurve].Tangent(u));
        last = curve[idCurve](u);
      }
    }
    else
    {
      break;
    }
  }

  // Append last point of the curve
  if (last != curve[curve.size() - 1](1.0))
  {
    p.push_back(curve[curve.size() - 1](1.0));  // Condition d'arret
    tangents.push_back(curve[curve.size() - 1].Tangent(0.9999));
  }
  return p;
}

/*!
\brief Generates a discretization of the piecewise cubic curve with a linear curvilign absisca parameterization.

\param step Stepping distance.
*/
std::vector<Vector2> CubicCurve2Set::GetDiscretisation(const double& step) const
{
  double niveau = 0.0;
  int idCurve = 0;

  // Set of points
  std::vector<Vector2> p;
  Vector2 last;

  // Start with the very first vertex of the piecewise curve
  p.push_back(curve[0](0));
  last = curve[0](0);

  //int parcours = 0;                                 // Condition d'arret
  while (true)
  {
    niveau += step;                              // Position suivante sur la trajectoire

    double u;
    idCurve = U(niveau, u);

    if (idCurve != -1)
    {
      if (last != curve[idCurve](u))
      {
        p.push_back(curve[idCurve](u));  // Condition d'arret
        last = curve[idCurve](u);
      }
    }
    else
    {
      break;
    }
  }

  // Append last point of the curve
  if (last != curve[curve.size() - 1](1.0))
  {
    p.push_back(curve[curve.size() - 1](1.0));  // Condition d'arret
  }
  return p;
}

/*!
\brief Compute the parameter of the curve corresponding to the input length.
\param s Input length.
\param u Parameter of the i-th curve.
\return Identifier of the curve in the piecewise definition.
*/
int CubicCurve2Set::U(const double& s, double& u) const
{
  if ((s > length) || (s < 0.0))
  {
    u = 0.0;
    return -1;
  }

  // Cumulative length
  double l = 0.0;

  for (int i = 0; i < curve.size(); i++)
  {
    // We found the right interval
    if (s < l + lengths[i])
    {
      // Find the right location
      u = curve[i].U(s - l, 256);
      return i;
    }
    // Increment the length
    l += lengths[i];
  }

  // Extreme case : we reached the end of the curve
  u = 1.0;
  return int(curve.size()) - 1;
}
