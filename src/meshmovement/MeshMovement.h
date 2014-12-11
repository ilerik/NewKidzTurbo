#ifndef TURBO_MeshMovement_MeshMovement
#define TURBO_MeshMovement_MeshMovement

#include "grid.h"
#include "basetypes.h"
#include "meshquality.h"
#include <cassert>

//Impelementation of mesh movement algorithms
class MeshMovement {
public:
	void IDWMove(Grid& grid, std::vector<int> nodes, std::vector<Vector> displacements) {
		assert(displacements.size() == nodes.size());

		for (int i = 0; i<nodes.size(); i++) {
			int nodeIndex = nodes[i];
				
		};


	};
};

#endif