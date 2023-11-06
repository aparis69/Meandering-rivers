#include "meanders.h"

// Simulation constants
double MeanderSimulation::Omega = -1.0;
double MeanderSimulation::Gamma = 2.5;
double MeanderSimulation::K = 1.0;
double MeanderSimulation::Cf = 0.011;
double MeanderSimulation::K1 = 60.0 / (365.0 * 24.0 * 60.0 * 60.0);
double MeanderSimulation::Dt = 9460800.0;
double MeanderSimulation::MaxSlope = 0.1;
double MeanderSimulation::Kv = 1.0e-12;
double MeanderSimulation::Dens = 1000;
double MeanderSimulation::tAvulsion = 1e-6;
double MeanderSimulation::pAvulsion = 0.01;
double MeanderSimulation::ChannelFalloff = 0.25;
double MeanderSimulation::SamplingDistance = 50.0;

/*!
\brief
*/
MeanderSimulation::MeanderSimulation()
{
	terrain = ScalarField2D(256, 256, Box2D(Vector2(0.0), 2500.0));
}

/*!
\brief
*/
MeanderSimulation::MeanderSimulation(const ScalarField2D& hf)
	: terrain(hf)
{
}

/*!
\brief
*/
Box2D MeanderSimulation::GetBox() const
{
	return terrain.GetBox();
}

/*!
\brief
*/
void MeanderSimulation::AddChannel(const Channel& ch)
{
	channels.push_back(ch);
}

/*!
\brief
*/
void MeanderSimulation::Step()
{
	ComputeMigrationRates();
	MigrateAllChannels();
	ManageCutoffs();
	ManageAvulsion();
	ResampleChannels();
}

/*!
\brief
*/
void MeanderSimulation::Step(int n)
{
	for (int i = 0; i < n; i++)
		Step();
}

/*!
\brief
*/
void MeanderSimulation::ComputeMigrationRates()
{
	for (int i = 0; i < channels.size(); i++)
		channels[i].ComputeMigrationRates();
}

/*!
\brief
*/
void MeanderSimulation::MigrateAllChannels()
{
	std::vector<Channel> results = channels;
	for (int i = 0; i < results.size(); i++)
	{
		if (results[i].Size() > 3)
			results[i].Migrate(GetBox(), terrain);
	}
	channels = results;
}

/*!
\brief
*/
void MeanderSimulation::ManageCutoffs()
{
	for (int i = 0; i < channels.size(); i++)
	{
		if (channels[i].Size() <= 3)
			continue;
		std::vector<Vector2>& points = channels[i].Points();
		for (int j = 0; j < points.size() - 2; j++)
		{
			Vector2 a = points[j];
			if (!terrain.Inside(a))
				continue;
			int n = int(points.size());
			int cutoffIndex = -1;
			for (int k = j + 5; k < n; k++)
			{
				Vector2 b = points[k];
				if (!terrain.Inside(b))
					continue;
				if (Magnitude(a - b) < channels[i].Width())
				{
					cutoffIndex = k;
					break;
				}
			}

			if (cutoffIndex != -1)
			{
				// Store oxbow lake for future visualization
				std::vector<Vector2> oxbowPts = channels[i].DoCutoff(cutoffIndex, j);

				// Advance main loop
				j = j + 1;
			}
		}
	}
}

/*!
\brief
*/
void MeanderSimulation::ManageAvulsion()
{
	for (auto& sec : channels)
	{
		// Check that the section is not too small
		// 50m x 50 points = 2.5km minimum for an avulsion to occur
		if (sec.Size() < 50)
			continue;

		const int padding = 10;
		for (int i = 1; i < sec.Size() - padding; i++)
		{
			// Avulsion probability is a function of total migration rate
			if (sec.MigrationRate(i) < tAvulsion)
				continue;
			if (Random::Uniform() < pAvulsion)
			{
				// Generate the new path
				std::vector<Vector2> path = sec.DoAvulsion(i, terrain);

				// Carve terrain along generated path to ensure correct flow
				//EnsureCorrectFlow(path);

				// Only one avulsion per step on a given section
				break;
			}
		}
	}
}

/*!
\brief
*/
void MeanderSimulation::ResampleChannels()
{
	for (int i = 0; i < channels.size(); i++)
	{
		if (channels[i].Size() > 3)
			channels[i].Resample();
	}
}
