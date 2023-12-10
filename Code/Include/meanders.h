#pragma once

#include "basics.h"
#include "curve.h"

#include <string>
#include <chrono>

struct MyChrono
{
public:
	std::chrono::time_point<std::chrono::high_resolution_clock> chrono;

	inline MyChrono()
	{
		chrono = std::chrono::high_resolution_clock::now();
	}
	inline void Restart()
	{
		chrono = std::chrono::high_resolution_clock::now();
	}
	inline long long ElapsedSeconds()
	{
		return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - chrono).count();
	}
	inline long long ElapsedMs()
	{
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - chrono).count();
	}
};

class Channel
{
private:
	std::vector<Vector2> pts;
	std::vector<double> ptsLocalMigrationRates;
	std::vector<double> ptsMigrationRates;
	double width;
	double depth;

public:
	inline Channel() { }
	Channel(const std::vector<Vector2>& pts, double w, double d);

	inline std::vector<Vector2>& Points() { return pts; }
	inline const std::vector<Vector2>& Points() const { return pts; }
	inline Vector2 Point(int i) const { return pts[i]; }
	inline double MigrationRate(int i) const { return ptsMigrationRates[i]; }
	inline double Width() const { return width; }

	int Size() const;
	double Length() const;
	double CurvilinearLength() const;
	double Sinuosity();
	Vector2 Tangent(int i) const;
	Vector2 MigrationDirection(int i) const;
	double Curvature(int i) const;
	double ScaledCurvature(int i) const;
	CubicCurve2Set ToCubicCurve() const;

	void Resample();
	void ComputeMigrationRates();
	void Migrate(const Box2D& domain, const ScalarField2D& terrain);
	std::vector<Vector2> DoCutoff(int cutoffIndex, int j);
	std::vector<Vector2> DoAvulsion(int startIndex, const ScalarField2D& bedrock);

private:
	bool Intersect(const Segment2& s, int startIndex, Vector2& hit, int& hitIndex) const;
	void ComputeLocalMigrationRates();
	void ComputeTotalMigrationRates();
	std::vector<Vector2> GeneratePath(int startIndex, int endIndex) const;
};

class MeanderSimulation
{
private:
	ScalarField2D terrain;
	std::vector<Channel> channels;

public:
	// Simulation parameters
	static double Omega;					//!< Constant in migration rate calculation (Howard and Knutson, 1984)
	static double Gamma;					//!< Constant from Ikeda et al., 1981 and Howard and Knutson, 1984
	static double K;							//!< Constant in Howard 1984 equation
	static double K1;							//!< Migration rate constant (m/s)
	static double Cf;							//!< Dimensionless Chezy friction factor
	static double Dt;							//!< Delta time (s)
	static double Kv;							//!< Vertical slope-dependent erosion rate constant (m/s)
	static double MaxSlope;				//!< Maximum terrain slope
	static double tAvulsion;			//!< Avulsion migration rate threshold.
	static double tAvulsionLength;//!< Minimum channel size for an avulsion to occur.
	static double ChannelFalloff;	//!< Channel falloff for start and end parts, in [0, 1]
	static double SamplingDistance;	//!< Maximum distance between points in a channel, in meters.

public:
	MeanderSimulation();
	MeanderSimulation(int seed);
	MeanderSimulation(int seed, const ScalarField2D& hf);

	// User control
	void AddChannel(const Channel& ch);
	void TriggerAvulsion();

	// Simulation
	void Step();
	void Step(int n);

	// Utility
	Box2D GetBox() const;
	void OutputImage(const std::string& path, int width, int height) const;

private:
	void EnsureCoherentFlow(const std::vector<Vector2>& path);
	void ComputeMigrationRates();
	void MigrateAllChannels();
	void ManageCutoffs();
	void ManageAvulsion();
	void ResampleChannels();
	bool SanityCheckChannels(const char* checkedFunction);
};
