#include "meanders.h"
#include "grid2.h"
#include <assert.h>

// Simulation constants
double MeanderSimulation::Omega = -1.0;
double MeanderSimulation::Gamma = 2.5;
double MeanderSimulation::K = 1.0;
double MeanderSimulation::Cf = 0.011;
double MeanderSimulation::K1 = 60.0 / (365.0 * 24.0 * 60.0 * 60.0);
double MeanderSimulation::Dt = 9460800.0;
double MeanderSimulation::MaxSlope = 0.1;
double MeanderSimulation::Kv = 1.0e-12;
double MeanderSimulation::tAvulsion = 1e-6;
double MeanderSimulation::tAvulsionLength = 2500.0;
double MeanderSimulation::ChannelFalloff = 0.05;
double MeanderSimulation::SamplingDistance = 50.0;

/*!
\brief Default constructor.
*/
MeanderSimulation::MeanderSimulation()
{
	terrain = ScalarField2D(256, 256, Box2D(Vector2(0.0), 2500.0));
	srand(1234);
}

/*!
\brief Constructor from a random seed.
\param seed the random seed
*/
MeanderSimulation::MeanderSimulation(int seed)
{
	terrain = ScalarField2D(256, 256, Box2D(Vector2(0.0), 2500.0));
	srand(seed);
}

/*!
\brief Constructor from a random seed and a given terrain.
\param seed the random seed
\param hf the terrain
*/
MeanderSimulation::MeanderSimulation(int seed, const ScalarField2D& hf)
	: terrain(hf)
{
	srand(seed);
}

/*!
\brief Returns the terrain bounding box.
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

	// Ensure that points are not perfectly aligned
	auto& pts = channels[channels.size() - 1].Points();
	for (int i = 1; i < pts.size() - 1; i++)
	{
		double xx = Random::Uniform(-2.5, 2.5);
		double yy = Random::Uniform(-2.5, 2.5);
		pts[i].x += xx;
		pts[i].y += yy;
	}

	// Resample the channel
	channels[channels.size() - 1].Resample();
}

/*!
\brief Manually trigger an avulsion event within a random channel in the network.
The starting point is choosen randomly within the three points that have the highest
curvature in absolute value.
*/
void MeanderSimulation::TriggerAvulsion()
{
	// Store [channel index, point index in channel]
	std::vector<std::pair<int, int>> candidatePointsInSections;
	for (int i = 0; i < channels.size(); i++)
	{
		// Check that the section is not too small
		if (channels[i].Length() < tAvulsionLength)
			continue;
		const int padding = 10;
		for (int j = 1; j < channels[i].Size() - padding; j++)
		{
			if (Math::Abs(channels[i].MigrationRate(j)) < tAvulsion)
				continue;
			candidatePointsInSections.push_back({ i, j });
		}
	}

	if (candidatePointsInSections.empty())
		return;

	// Stochastically choose starting index in all candidates
	const int randomIndex = Random::Integer(int(candidatePointsInSections.size()));
	const std::pair<int, int> avulsionData = candidatePointsInSections[randomIndex];

	std::vector<Vector2> newPath = channels[avulsionData.first].DoAvulsion(avulsionData.second, terrain);
	EnsureCoherentFlow(newPath);
	ResampleChannels();
}

/*!
\brief
*/
void MeanderSimulation::Step()
{
	ComputeMigrationRates();
	assert(SanityCheckChannels("ComputeMigrationRates"));

	MigrateAllChannels();
	assert(SanityCheckChannels("MigrateAllChannels"));

	ManageCutoffs();
	assert(SanityCheckChannels("ManageCutoffs"));

	ResampleChannels();
	assert(SanityCheckChannels("ResampleChannels"));
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
\brief Ensure a coherent flow in the terrain along a given path.
Put simply, this function makes sure that the elevation is decreasing along the trajectory.
\param path the trajectory
*/
void MeanderSimulation::EnsureCoherentFlow(const std::vector<Vector2>& path)
{
	const double radius = Magnitude(terrain.CellDiagonal());
	for (int i = 1; i < path.size(); i++)
	{
		double d = terrain.GetValueBilinear(path[i - 1]) - terrain.GetValueBilinear(path[i]);
		if (d < 0)
		{
			int a, b;
			terrain.VertexToInteger(path[i], a, b);
			terrain.Set(a, b, terrain.GetValueBilinear(path[i - 1]) - 5.0);
		}
	}
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
void MeanderSimulation::ResampleChannels()
{
	for (int i = 0; i < channels.size(); i++)
	{
		if (channels[i].Size() > 3)
			channels[i].Resample();
	}
}

/*!
\brief
*/
void MeanderSimulation::OutputImage(const std::string& path, int width, int height) const
{
	auto rasterizeBoxFunc = [](Grid2<Vector3>& img, int px, int pxx, int py, int pyy)
	{
		int px0 = px;
		int py0 = py;
		for (px = px0; px < Math::Min(pxx, img.Width()); px++)
		{
			for (py = py0; py < Math::Min(pyy, img.Height()); py++)
				img(px, py) = Vector3(1.0, 1.0, 1.0);
		}
	};

	const Box2D domain = GetBox();
	const ScalarField2D dummy = ScalarField2D(width, height, domain);
	Grid2<Vector3> image(width, height, Vector3(0.0));
	for (const auto& ch : channels)
	{
		// Rasterize all points as small boxes
		const auto& pts = ch.Points();
		for (int i = 0; i < pts.size(); i++)
		{
			int px, py;
			dummy.VertexToInteger(pts[i], px, py);
			rasterizeBoxFunc(image, px, px + 2, py, py + 2);
		}
	}

	// Export to ppm file
	FILE* fp = NULL;
	fp = fopen(path.c_str(), "wb");
	if (fp == NULL)
	{
		std::cout << "Couldn't write to file - exiting" << std::endl;
		return;
	}
	fprintf(fp, "P6\n%d %d\n255\n", image.Width(), image.Height());
	for (int i = 0; i < image.Width(); i++)
	{
		for (int j = 0; j < image.Height(); j++)
		{
			static unsigned char color[3];
			Vector3 c = image(i, j);
			color[0] = ((int)c[0] * 255) % 256;
			color[1] = ((int)c[1] * 255) % 256;
			color[2] = ((int)c[2] * 255) % 256;
			(void)fwrite(color, 1, 3, fp);
		}
	}
	fclose(fp);
}

/*!
\brief
*/
bool MeanderSimulation::SanityCheckChannels(const char* checkedFunction)
{
	for (auto& channel : channels)
	{
		for (int i = 0; i < channel.Size(); i++)
		{
			if (Math::IsNumber(channel.Point(i)[0]) == false
				|| Math::IsNumber(channel.Point(i)[1]) == false
				|| Math::IsNumber(channel.MigrationRate(i)) == false)
			{
				std::cout << checkedFunction << ": NaN error" << std::endl;
				assert(false);
				return false;
			}
		}
	}
	return true;
}
