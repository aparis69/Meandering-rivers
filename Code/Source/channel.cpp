#include "meanders.h"

Channel::Channel(const std::vector<Vector2>& pts, double w, double d)
	: pts(pts), width(w), depth(d)
{

}

Vector2 Channel::Tangent(int i) const
{
	const int k = pts.size() - 1;
	if (i == 0)
		return pts.at(1) - pts.at(0);
	else if (i == k)
		return pts.at(k) - pts.at(k - 1);
	else
		return (pts.at(i + 1) - pts.at(i - 1)) / 2.0;
}

double Channel::Curvature(int i) const
{
	if (i == 0 || i == pts.size() - 1)
		return 0.0;

	// First derivative
	Vector2 dxy = Tangent(i);

	// Second derivative
	Vector2 dxy0 = Tangent(i - 1);
	Vector2 dxy1 = Tangent(i + 1);
	Vector2 ddxy = (dxy1 - dxy0) / 2.0;

	// Compute curvature
	double dx = dxy[0];
	double dy = dxy[1];
	double ddx = ddxy[0], ddy = ddxy[1];
	return (dx * ddy - dy * ddx) / Math::Pow(Math::Sqr(dx) + Math::Sqr(dy), 1.5);
}

double Channel::ScaledCurvature(int i) const
{
	return width * Curvature(i);
}

void Channel::Resample() const
{

}
