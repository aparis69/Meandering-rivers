#pragma once

#include "basics.h"

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
	double width;
	double depth;

public:
	inline Channel() { }
	Channel(const std::vector<Vector2>& pts, double w, double d);
	
	Vector2 Tangent(int i) const;
	double Curvature(int i) const;
	double ScaledCurvature(int i) const;
	void Resample() const;
};

class Network 
{
private:
	ScalarField2D terrain;
	std::vector<Channel> channels;

public:
	// Simulation parameters
	static double Omega;			//!< Constant in migration rate calculation (Howard and Knutson, 1984)
	static double Gamma;			//!< Constant from Ikeda et al., 1981 and Howard and Knutson, 1984
	static double K;				//!< Constant in Howard 1984 equation
	static double K1;				//!< Migration rate constant (m/s)
	static double Cf;				//!< Dimensionless Chezy friction factor
	static double Dt;				//!< Delta time (s)
	static double Kv;				//!< Vertical slope-dependent erosion rate constant (m/s)
	static double Dens;				//!< Density of water (kg/m3)
	static double MaxSlope;			//!< Maximum slope
	static double tAvulsion;		//!< Avulsion local curvature threshold.
	static double pAvulsion;		//!< Avulsion event probability (given that curvature > tAvulsion).
	static double ChannelFalloff;	//!< Channel falloff for start and end parts, in [0, 1]
	static double SamplingDistance;	//!< Maximum distance between points in a channel, in meters.

public:
	Network();
	Network(const ScalarField2D& hf);
	
	void AddChannel(const Channel& ch);
	void Step();
	void Step(int n);

private:
	void ComputeMigrationRates();
	void MigrateChannels();
	void ManageCutoffs();
	void ManageAvulsion();
	void ResampleChannels();
};
