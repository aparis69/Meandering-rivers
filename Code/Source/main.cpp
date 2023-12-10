/*
 * This is a reimplementation of the paper "Authoring and simulating meandering rivers" published at Transactions on graphics
 * and presented at Siggraph Asia 2023. Running this will output a serie of .ppm files showing the simulation output (a set of curves, basically)
 *
 * The code has minimal dependencies and no real time visualization is provided, so it should be straightforward to use, compile, and run.
 *
 * If you have any question or problem to compile the code, you can contact me at:
 * axel.paris69@gmail.com
*/

#include "meanders.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

/*!
\brief Perform the simulation for a given number of steps, and optionally rasterize and export the simulation frames.
At the very least, the last frame is saved to the disk as a PPM file.
\param simu the simulation, must be properly initialized
\param N step count
\param exportAllSteps if true, will export all simulation steps in PPM files
\param w, h dimensions of the rasterized frames
*/
static void Simulate(MeanderSimulation& simu, int N, bool exportAllSteps = false, int w = 1024, int h = 1024)
{
	if (exportAllSteps)
		simu.OutputImage("../Results/step0.ppm", w, h);
	MyChrono timer;
	for (int i = 1; i <= N; i++)
	{
		simu.Step(1);
		if (exportAllSteps)
			simu.OutputImage("../Results/step" + std::to_string(i) + ".ppm", w, h);
	}
	std::cout << "Performed " << N << " simulation steps in " << timer.ElapsedMs() << "ms" << std::endl;
	if (!exportAllSteps)
		simu.OutputImage("../Results/step" + std::to_string(N) + ".ppm", w, h);
}

/*!
\brief Initialize a heightfield from a png file.
\param terrain the resulting heightfield
\param filename relative file path
*/
static ScalarField2D LoadScalarFieldFromFile(double cellSizeX, double cellSizeY, double maxZ, const char* filename)
{
	int x, y, n;
	unsigned char* data = stbi_load(filename, &x, &y, &n, 0);
	ScalarField2D terrain = ScalarField2D(x, y, Box2D(x * cellSizeX, y * cellSizeY));
	for (int i = 0; i < x * y; i += n)
	{
		terrain.Set(i, ((double)data[i]) * maxZ);
	}
	stbi_image_free(data);
	return terrain;
}

/*!
\brief Perform a basic simulation over a single channel.
Output is one PPM file showing the final state of the channel.
*/
static void ExampleSingleChannel()
{
	// Init
	ScalarField2D terrain(256, 256, Box2D(Vector2(0), 10000.0));
	std::vector<Vector2> pts =
	{
		Vector2(-7500, 0),
		Vector2(-1750, 0),
		Vector2(-1000, 0),
		Vector2(0, 0),
		Vector2(1000, 0),
		Vector2(1500, 0),
		Vector2(7500, 0)
	};
	MeanderSimulation simulation(1234, terrain);
	simulation.AddChannel(Channel(pts, 50.0, 2.0));

	// Simulate
	Simulate(simulation, 500);
}

/*!
\brief Perform some simulation step, and then trigger an avulsion event manually.
Output is two PPM files.
*/
static void ExampleAvulsion()
{
	// Init
	ScalarField2D terrain(256, 256, Box2D(Vector2(0), 10000.0));
	std::vector<Vector2> pts =
	{
		Vector2(-7500, 0),
		Vector2(-1750, 0),
		Vector2(-1000, 0),
		Vector2(0, 0),
		Vector2(1000, 0),
		Vector2(1500, 0),
		Vector2(7500, 0)
	};
	MeanderSimulation simulation(1234, terrain);
	simulation.AddChannel(Channel(pts, 50.0, 2.0));

	// Simulate
	simulation.Step(400);
	simulation.OutputImage("../Results/avulsion_before.ppm", 1024, 1024);
	simulation.TriggerAvulsion();
	simulation.Step(1);
	simulation.OutputImage("../Results/avulsion_after.ppm", 1024, 1024);
}

/*!
\brief Perform a simulation constrained by a given terrain file (loaded from disk).
Output is one PPM file showing the final state of the channel. The terrain is also shaded
to see its influence.
*/
static void ExampleTerrainConstrained()
{
	// Init
	ScalarField2D terrain = LoadScalarFieldFromFile(125.0, 137.0, 2500.0, "../Resources/hf.png");

	std::vector<Vector2> pts =
	{
		Vector2(-6400.46, -16602.8),
		Vector2(-5778.78, -15109.9),
		Vector2(-5256.07, -13642.9),
		Vector2(-4437.77, -11791.4),
		Vector2(-4153.77, -10096.8),
		Vector2(-3791.16, -8568.63),
		Vector2(-3248.14, -7130.27),
		Vector2(-2654.12, -5671.92),
		Vector2(-2071.26, -3966.19),
		Vector2(-1761.29, -1759.61),
		Vector2(-1496.86, -66.0183),
		Vector2(-1571.23, 1797.13),
		Vector2(-1708.55, 3404.22),
		Vector2(-1891.34, 4656.92),
		Vector2(-2386.36, 5772.93),
		Vector2(-3053.14, 6754.97),
		Vector2(-4101.77, 7343.27),
		Vector2(-5098.04, 7681.74),
		Vector2(-6419.52, 7983.41),
		Vector2(-7571.98, 8250.98),
		Vector2(-8649.53, 8283.52),
		Vector2(-9524.27, 8273.19)
	};
	MeanderSimulation simulation(1234, terrain);
	simulation.AddChannel(Channel(pts, 100.0, 4.0)); // TODO: check automatic depth computation

	// Simulate
	simulation.Step(400);

	// TODO: fix output image
	simulation.OutputImage("../Results/meander_constrained.ppm", terrain.SizeY() * 10, terrain.SizeX() * 10);
}

/*!
\brief Perform a more complex simulation involving multiple channels
as well as an underlying terrain. Initial river network was computed offline,
and is a mix of "automatically computed from the heightfield" and "manually set by a user".
*/
static void ExampleNetwork()
{
	// TODO
}

int main()
{
	ExampleSingleChannel();
	ExampleAvulsion();

	// These examples are not yet complete/working :-)
	//ExampleTerrainConstrained();
	//ExampleNetwork();

	return 0;
}
