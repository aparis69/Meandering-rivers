#include "meanders.h"

Channel::Channel(const std::vector<Vector2>& pts, double w)
	: pts(pts), width(w)
{
	// River depth can be deduced from width (Konsoer et al. 2013, see Sylvester 2019 supplementary material)
	depth = pow(width / 18.8, 1.0 / 1.41);
	depth = Math::Max(1.0, depth);
	ptsMigrationRates.resize(pts.size(), 0.0);
}

int Channel::Size() const
{
	return int(pts.size());
}

double Channel::Length() const
{
	return Magnitude(pts[0] - pts[pts.size() - 1]);
}

double Channel::CurvilinearLength() const
{
	double length = 0.0;
	for (int i = 0; i < pts.size() - 1; i++)
	{
		length += Magnitude(pts.at(i + 1) - pts.at(i));
	}
	return length;
}

double Channel::Sinuosity()
{
	return CurvilinearLength() / Length();
}

Vector2 Channel::Tangent(int i) const
{
	const int k = int(pts.size()) - 1;
	if (i == 0)
		return pts.at(1) - pts.at(0);
	else if (i == k)
		return pts.at(k) - pts.at(k - 1);
	else
		return (pts.at(i + 1) - pts.at(i - 1)) / 2.0;
}

Vector2 Channel::Normal(int i) const
{
	return Tangent(i).Orthogonal();
}

Vector2 Channel::MigrationDirection(int i, const std::vector<PointConstraint>& constraints) const
{
	// Base direction is normal to the curve
	Vector2 dir = -Normalize(Normal(i));

	// Modify based on user constraints
	const Vector2 p = Point(i);
	double nb = 1.0;
	for (int k = 0; k < constraints.size(); k++)
	{
		Vector2 d = constraints[k].Direction(p);
		if (d != Vector2(0, 0))
		{
			dir += d;
			nb += 1.0;
		}
	}
	dir /= nb;
	return dir;
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

CubicCurve2Set Channel::ToCubicCurve() const
{
	return CubicCurve2Set(pts);
}


void Channel::Resample()
{
	// Resample along a cubic curve
	const double x = MeanderSimulation::SamplingDistance;
	CubicCurve2Set cubicCurve = ToCubicCurve();
	std::vector<Vector2> newPts = cubicCurve.GetDiscretisation(x);

	// Copy to existing points
	const int count = int(Math::Min(newPts.size(), pts.size() - 1));
	for (int j = 1; j < count; j++)
	{
		double z = pts[j][2];
		pts[j] = newPts[j];
	}

	// Add new points if necessary
	for (int j = int(pts.size()) - 1; j < int(newPts.size()) - 1; j++)
	{
		pts.insert(pts.begin() + j, newPts[j]);
	}

	// Or remove last ones
	const int pointCountToRemove = int(pts.size() - newPts.size());
	for (int j = 0; j < pointCountToRemove; j++)
		pts.erase(pts.end() - 1);

	ptsLocalMigrationRates.resize(pts.size(), 0.0);
	ptsMigrationRates.resize(pts.size(), 0.0);
}

void Channel::ComputeMigrationRates()
{
	ComputeLocalMigrationRates();
	ComputeTotalMigrationRates();
}

void Channel::Migrate(
	const Box2D& domain, 
	const Grid2<Vector2>& terrainGrad,
	const std::vector<PointConstraint>& constraints)
{
	std::vector<Vector2> outPoints = pts;
	const double channelSize = Math::Sqr(Length());
	const double MaxSlopeSqr = Math::Sqr(MeanderSimulation::MaxSlope);
	for (int i = 1; i < pts.size() - 1; i++)
	{
		if (!domain.Contains(pts[i]))
			continue;

		// Falloff at start/end of section
		double d1 = SquaredMagnitude(pts[i] - pts[0]);
		double d2 = SquaredMagnitude(pts[i] - pts[pts.size() - 1]);
		double wf = 1.0 - Math::CubicSmoothCompact(Math::Min(d1, d2), channelSize * MeanderSimulation::ChannelFalloff);

		// Falloff with terrain slope
		Vector2 q = pts[i] - domain[0];
		double u = q[0] / (domain[1].x - domain[0].x);
		double v = q[1] / (domain[1].y - domain[0].y);
		double sq = SquaredMagnitude(terrainGrad.BilinearValue(v, u));
		double wt = Math::CubicSmoothCompact(sq, MaxSlopeSqr);

		// Combine all falloffs
		double w = wf * wt;

		// Point migration
		Vector2 normalMigration = MigrationDirection(i, constraints);
		outPoints[i][0] = pts[i][0] + w * MeanderSimulation::Dt * ptsMigrationRates[i] * normalMigration[0];
		outPoints[i][1] = pts[i][1] + w * MeanderSimulation::Dt * ptsMigrationRates[i] * normalMigration[1];
	}
	pts = outPoints;
}

std::vector<Vector2> Channel::DoCutoff(int cutoffIndex, int j)
{
	// points[j] must be linked to points[cutoffIndex]
	// Delete everything inbetween
	std::vector<Vector2> oxbowPts;
	int limit = cutoffIndex;
	int deleteIndex = j + 1;
	oxbowPts.push_back(pts[j]);
	for (int k = j + 1; k < limit; k++)
	{
		oxbowPts.push_back(pts[deleteIndex]);
		pts.erase(pts.begin() + deleteIndex);
		cutoffIndex--;
	}
	oxbowPts.push_back(pts[cutoffIndex]);
	return oxbowPts;
}

std::vector<Vector2> Channel::DoAvulsion(int startIndex, const ScalarField2D& bedrock)
{
	// Create new path between start and end
	int endIndex = startIndex + 1 + Random::Integer(int(pts.size() - startIndex - 1));
	std::vector<Vector2> avulsionPath = GeneratePath(startIndex, endIndex);

	// Create new vector of points for the section
	std::vector<Vector2> newPoints;
	{
		// (1) Add points before the avulsion, which are not modified
		for (int i = 0; i < startIndex; i++)
			newPoints.push_back(pts[i]);

		// (2) Add points for the avulsion, and check intersection with existing path
		for (int i = 0; i < avulsionPath.size(); i++)
		{
			if (i < avulsionPath.size() - 1)
			{
				Segment2 part = Segment2(avulsionPath[i], avulsionPath[i + 1]);
				Vector2 hit;
				int hitIndex;
				if (Intersect(part, startIndex + 1 + i, hit, hitIndex))
				{
					avulsionPath.erase(avulsionPath.begin() + i, avulsionPath.end());
					avulsionPath.push_back(hit);

					newPoints.insert(newPoints.begin() + startIndex + i, avulsionPath[i]);

					endIndex = hitIndex;
					break;
				}
			}

			newPoints.insert(newPoints.begin() + startIndex + i, avulsionPath[i]);
		}

		// (3) Add points after the avulsion, which are not modified
		for (int i = endIndex; i < pts.size(); i++)
			newPoints.push_back(pts[i]);
	}

	// Assign!
	pts = newPoints;

	// Return avulsion path
	return avulsionPath;
}


bool Channel::Intersect(const Segment2& s, int startIndex, Vector2& hit, int& hitIndex) const
{
	for (int i = startIndex; i < pts.size() - 1; i++)
	{
		Segment2 part = Segment2(pts[i], pts[i + 1]);
		if (part.Intersect(s, hit))
		{
			hitIndex = i + 1;
			return true;
		}
	}
	hit = Vector2(0, 0);
	return false;
}

void Channel::ComputeLocalMigrationRates()
{
	for (int i = 0; i < pts.size(); i++)
		ptsLocalMigrationRates[i] = MeanderSimulation::K1 * ScaledCurvature(i);
}

void Channel::ComputeTotalMigrationRates()
{
	if (pts.size() <= 3)
		return;

	// Precompute river section lengths
	std::vector<double> lengths;
	lengths.resize(pts.size(), 0.0);
	double l = 0.0;
	for (int i = 0; i < lengths.size() - 1; i++)
	{
		lengths[i] = l;
		l = l + Magnitude(pts[i] - pts[i + 1]);;
	}
	lengths[lengths.size() - 1] = l;

	// Final migration rate (function of local migration rate and sinuosity)
	const double sinuosity = Math::Pow(Sinuosity(), -2.0 / 3.0);
	const double alpha = MeanderSimulation::K * 2.0 * MeanderSimulation::Cf / depth;
	for (int i = 0; i < pts.size(); i++)
	{
		// This value drastically change the output of the simulation.
		// In geomorphology, we go toward the starting point of the river
		// And compute a contribution (minimal, but not zero).
		// This version is faster but gives different results (speedup of ~25%)
		int start = i - 1500;
		start = Math::Clamp(start, 0, start);
		int size = i - start;

		std::vector<double> si2;
		std::vector<double> G;
		si2.resize(size + 1);
		G.resize(size + 1);
		for (int k = 0; k < G.size(); k++)
		{
			si2[G.size() - 1 - k] = lengths[i - k];
			G[G.size() - 1 - k] = exp(-alpha * si2[G.size() - 1 - k]);
		}

		double sumR0 = 0.0;
		double sumG = 0.0;
		for (int k = 0; k < G.size(); k++)
		{
			sumR0 += ptsLocalMigrationRates[i - k] * G[k];
			sumG += G[k];
		}

		// Avoid NaN issues!
		if (sumG == 0.0)
			sumG = 1.0;

		ptsMigrationRates[i] = MeanderSimulation::Omega * ptsLocalMigrationRates[i] + MeanderSimulation::Gamma * sumR0 / sumG;
		ptsMigrationRates[i] = sinuosity * ptsMigrationRates[i];
	}
}

std::vector<Vector2> Channel::GeneratePath(int startIndex, int endIndex) const
{
	const Vector2 start = pts[startIndex];
	const Vector2 end = pts[endIndex];

	Vector2 dir = Normalize(end - start);
	double angle = Random::Uniform(-35.0, 35.0);
	Vector2 perturbDir = RotationMatrix(Math::DegreeToRadian(angle)) * dir;

	Vector2 p = start;
	const double totalLength = SquaredMagnitude(start - end);
	std::vector<Vector2> path;
	while (SquaredMagnitude(p - end) > Math::Sqr(MeanderSimulation::SamplingDistance))
	{
		dir = Normalize(end - p);
		double t = 1.0 - (SquaredMagnitude(p - end) / totalLength);
		Vector2 d = Lerp(perturbDir, dir, t);
		p = p + d * MeanderSimulation::SamplingDistance;
		p = p + Vector2(Random::Uniform(-10.0, 10.0), Random::Uniform(-10.0, 10.0));
		path.push_back(p);
	}

	return path;
}
