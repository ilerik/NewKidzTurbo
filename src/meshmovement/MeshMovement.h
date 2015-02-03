#ifndef TURBO_MeshMovement_MeshMovement
#define TURBO_MeshMovement_MeshMovement

#include "grid.h"
#include "basetypes.h"
#include "meshquality.h"
#include "geomfunctions.h"
#include <cassert>
#include <unordered_set>

//Impelementation of mesh movement algorithms
class MeshMovement {
public:
	//Weighting Function
	double W(Vector d_r) {	
		double dr = d_r.mod();
		double drx = d_r.x;
		if (abs(dr) < 1e-15) return 1;
		double Ai = 1;
		double Ldef = 1e-2;
		double alpha = 0.1;
		double a = 3.0;
		dr = d_r.mod();
		//double a = 2.5;
		double b = 5;
		double res = Ai *( pow(Ldef / dr, a) + pow(alpha * Ldef / dr, b));
		return res;
	};


	void IDWComputeDisplacements(Grid& grid, 
		const std::unordered_set<int> movingNodes, 
		std::map<int, Vector>& movingNodesDisplacements, 
		const std::unordered_set<int> freeNodes, 
		std::map<int, Vector>& freeNodesDisplacements) 
	{
		double lambdaX = 0.8 * 1e-2; //Wave number [cm]
		int ModesNumber = 1; //Number of modes
		double xMax = ModesNumber * (lambdaX * 0.5);
		double xMin = ModesNumber * (-lambdaX * 0.5);
		double period = xMax - xMin;

		//assert(movingNodesDisplacements.size() == movingNodes.size());
		//For all nodes allocate memory to store rotation and displacements
		std::map<int, Vector> dR;
		std::map<int, Vector> b;
		std::map<int, RotationMatrix> R;
		std::map<int, std::set<int> > bFaces;

		//For all nodes
		for (int movingNodeIndex : movingNodes) {
			if (movingNodesDisplacements.find(movingNodeIndex) == movingNodesDisplacements.end()) throw new Exception("Displacement for node is absent");
			dR[movingNodeIndex] = movingNodesDisplacements[movingNodeIndex];
		};
				
		//For each moving node determine neighbour boundary faces	
		for (int i = 0; i<grid.localFaces.size(); i++) {
			Face& f = grid.localFaces[i];
			bool isMovingFace = true;
			for (int& faceNodeIndex : f.FaceNodes) {
				if (movingNodes.find(faceNodeIndex) == movingNodes.end()) {
					isMovingFace = false;
					break;
				};
			};

			if (isMovingFace) for (int& faceNodeIndex : f.FaceNodes) {
				bFaces[faceNodeIndex].insert(i);
			};
		};	
		//std::cout<<"N="<<nodes.size()<<", Nb="<<bNodes.size()<<'\n';
	
		//Iterate through all moving nodes compute rotation
		int counter = 0;
		for (int nodeIndex : movingNodes) {
			counter++;		
			Node& node = grid.localNodes[nodeIndex];
			Vector Displacement = dR[nodeIndex];
			Vector newNormal;
			Vector normal;
			int k = 0;
			//Iterate through all neighbour faces
			std::set<int> faces = bFaces[nodeIndex];
			for (int faceIndex : faces) {
				Face& f = grid.localFaces[faceIndex];
				std::vector<Vector> points(0);
				for (int i = 0; i<f.FaceNodes.size(); i++) {
					Node& n = grid.localNodes[f.FaceNodes[i]];
					Vector D = dR[f.FaceNodes[i]];
					points.push_back(n.P + D);
				};
				Vector n = CalcNormal(points);
				if (f.FaceNormal * n > 0) {
					newNormal += n;
				} else {
					newNormal -= n;
				};			
				normal += f.FaceNormal / f.FaceNormal.mod();
				k++;
			};				
			Vector r = newNormal;
			newNormal /= newNormal.mod();		
			normal /= normal.mod();
	//		std::cout<<"Rotation : "<<counter<<" normal" << normal.mod() << "\n";
			RotationMatrix M = CalcRotation(normal, newNormal);
			R[nodeIndex] = M;
			b[nodeIndex] = node.P + Displacement - M * node.P;			
		};

		//Compute displacement for free nodes
		counter = 0;
		double percent = 0;
		for (int nodeIndex : freeNodes) {
			counter++;		
			Node& ni = grid.localNodes[nodeIndex];
			std::vector<Vector> displacements;
			dR[nodeIndex] = Vector(0,0,0);
			double sum = 0;		

			for (int movingNodeIndex : movingNodes) {			
				Node& nb = grid.localNodes[movingNodeIndex];
				Vector bb = b[movingNodeIndex];
				RotationMatrix M = R[movingNodeIndex];			
				Vector dr = ni.P - nb.P;			
				//Vector displ = (M * ni.P + bb - ni.P);
				Vector displ = dR[movingNodeIndex];
				/*if ((dr.x) > (period / 2.0)) {
					dr.x = (period / 2.0) - dr.x;
				};
				if ((displ.x) > (period / 2.0)) {
					displ.x = (period / 2.0) - displ.x;
				};*/
				double w = W(dr);	//TO DO		
				dR[nodeIndex] += w * displ; 
				sum += w; //			
			}					
			dR[nodeIndex] /= sum;	

			//Save result
			freeNodesDisplacements[nodeIndex] = dR[nodeIndex];
		};

		//Sync
	};

	void IDWMoveNodes(Grid& grid, std::vector<int> nodes, std::vector<Vector> displacements) {
		assert(displacements.size() == nodes.size());

		//Determine set of boundary nodes
		std::unordered_set<int> movingNodes(nodes.begin(), nodes.end());

		//Determine set of inner nodes		
		std::unordered_set<int> freeNodes;
		int cb = 0;
		for (int i = 0; i < grid.localNodes.size(); i++) {
			if (i == nodes[cb]) {
				cb++;
				continue;
			};
			freeNodes.insert(i);
		};		

		//For all nodes allocate memory to store rotation and displacements
		std::map<int, Vector> dRmoving;
		std::map<int, Vector> dRfree;

		//For boundary nodes
		for (int i = 0; i<nodes.size(); i++) {
			int movingNodeIndex = nodes[i];
			Vector displacement = displacements[i];
			dRmoving[movingNodeIndex] = displacement;
		};

		//Compute displacements
		IDWComputeDisplacements(grid, movingNodes, dRmoving, freeNodes, dRfree);

		//Change grid according to displacements	
		//Nodes first
		for (int i : freeNodes) {
			grid.localNodes[i].P += dRfree[i];
		};
		for (int i : movingNodes) {
			grid.localNodes[i].P += dRmoving[i];
		};

		//Recalculate Face normals and other grid geometric properties
		grid.UpdateGeometricProperties();
	};
};

#endif