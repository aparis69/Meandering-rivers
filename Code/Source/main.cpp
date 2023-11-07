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
	simulation.Step(350);
	simulation.OutputImage("../Results/avulsion_before.ppm", 1024, 1024);
	simulation.TriggerAvulsion();
	simulation.OutputImage("../Results/avulsion_after.ppm", 1024, 1024);
}

/*!
\brief Perform a simulation constrained by a given terrain file (loaded from disk).
Output is one PPM file showing the final state of the channel. The terrain is also shaded
to see its influence.
*/
static void ExampleTerrainConstrained()
{
	// TODO
}

/*!
\brief Perform a more complex simulation involving multiple channels
as well as an underlying terrain. Initial rivers were computed offline, 
and are a mix of automatically computed from the heightfield/manually set by a user.
*/
static void ExampleNetwork()
{
	// TODO
}

int main()
{
	ExampleSingleChannel();

	// These examples are not yet complete :-)
	//ExampleAvulsion();
	//ExampleTerrainConstrained();
	//ExampleNetwork();

	return 0;
}
