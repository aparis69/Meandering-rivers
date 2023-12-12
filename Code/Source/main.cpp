/*
 * This is a reimplementation of the paper "Authoring and simulating meandering rivers" published at Transactions on graphics
 * and presented at Siggraph Asia 2023. Running this will output a serie of .ppm files showing the simulation output (a set of curves, basically)
 *
 * The code has minimal dependencies and no real time visualization is provided, so it should be straightforward to use, compile, and run.
 * 
 * Keep in mind that these examples are kept simple on purposes, but more complex examples were made for the paper using the same algorithms and techniques.
 *
 * If you have any question or problem to compile the code, you can contact me at:
 * axel.paris69@gmail.com
*/

#include "meanders.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

/*!
\brief Initialize a heightfield from a png file.
\param terrain the resulting heightfield
\param filename relative file path
*/
static ScalarField2D LoadScalarFieldFromFile(double cellSizeX, double cellSizeY, double maxZ, const char* filename)
{
	int x, y, n;
	unsigned char* data = stbi_load(filename, &x, &y, &n, 1);
	ScalarField2D terrain = ScalarField2D(x, y, Box2D(x * cellSizeX, y * cellSizeY));
	int index = 0;
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < y; j++)
		{
			terrain.Set(i, j, (((double)data[index]) / 255.0) * maxZ);
			index++;
		}
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
	simulation.AddChannel(Channel(pts, 50.0));
	simulation.Step(500);
	simulation.OutputImage("../Results/meander_simple.ppm", 1024, 1024);
}

/*!
\brief Perform some simulation step, and then trigger an avulsion event manually.
Output is two PPM files.
*/
static void ExampleAvulsion()
{
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
	simulation.AddChannel(Channel(pts, 50.0));
	simulation.Step(700);
	simulation.OutputImage("../Results/avulsion_before.ppm", 1024, 1024);
	simulation.TriggerAvulsion();
	simulation.TriggerAvulsion();
	simulation.Step(50);
	simulation.OutputImage("../Results/avulsion_after.ppm", 1024, 1024);
}

/*!
\brief Perform a simulation constrained by a given terrain file (loaded from disk).
Output is one PPM file showing the final state of the channel. The terrain is also shaded
to see its influence.
*/
static void ExampleTerrainConstrained()
{
	ScalarField2D terrain = LoadScalarFieldFromFile(200.0, 200.0, 2500.0, "../Resources/hf.png");
	std::vector<Vector2> pts =
	{
		Vector2(0, -8500),
		Vector2(0, -1750),
		Vector2(0, -1000),
		Vector2(0, 0),
		Vector2(0, 1000),
		Vector2(0, 1500),
		Vector2(0, 8500)
	};
	MeanderSimulation simulation(1234, terrain);
	simulation.AddChannel(Channel(pts, 50.0));
	simulation.Step(1200);
	simulation.OutputImage("../Results/meander_terrain.ppm", 1024, 1024);
}

/*!
\brief Example of a single channel with an attractive point constraint.
*/
static void ExampleAttractiveConstraint()
{
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
	simulation.AddChannel(Channel(pts, 50.0));
	simulation.AddPointConstraint(PointConstraint(Vector2(0, 1500), 10000, 0.1));
	simulation.Step(500);
	simulation.OutputImage("../Results/meander_simple_constrained.ppm", 1024, 1024);
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
	//ExampleSingleChannel();
	ExampleAvulsion();
	ExampleTerrainConstrained();

	// These examples are not yet complete/working :-)
	//ExampleAttractiveConstraint();
	//ExampleNetwork();

	return 0;
}
