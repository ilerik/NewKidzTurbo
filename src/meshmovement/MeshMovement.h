#ifndef TURBO_MeshMovement_MeshMovement
#define TURBO_MeshMovement_MeshMovement

#include "grid.h"
#include "basetypes.h"
#include "meshquality.h"
#include "geomfunctions.h"
#include <cassert>

//Impelementation of mesh movement algorithms
class MeshMovement {
public:
	void IDWMove(Grid& grid, std::vector<int> nodes, std::vector<Vector> displacements) {
		assert(displacements.size() == nodes.size());

		//Determine set of boundary nodes
		std::set<int> bNodes(nodes.begin(), nodes.end());

		//Determine set of inner nodes		
		std::set<int> iNodes;
		int cb = 0;
		for (int i = 0; i < grid.localNodes.size(); i++) {
			if (i == nodes[cb]) {
				cb++;
				continue;
			};
			iNodes.insert(i);
		};		

		//For all nodes allocate memory to store rotation and displacements
		std::map<int, Vector> dR;
		std::map<int, Vector> b;
		std::map<int, RotationMatrix> R;
		std::map<int, std::set<int> > bFaces;

		//For all nodes
		for (int i = 0; i<nodes.size(); i++) {
			int nodeIndex = nodes[i];
			dR[nodeIndex] = displacements[i];
		};
				
		//For each moving node determine neighbour boundary faces	
		for (int i = 0; i<grid.localFaces.size(); i++) {
			Face& f = grid.localFaces[i];
			bool isMovingFace = true;
			for (int& faceNodeIndex : f.FaceNodes) {
				if (bNodes.find(faceNodeIndex) == bNodes.end()) {
					isMovingFace = false;
					break;
				};
			};

			if (isMovingFace) for (int& faceNodeIndex : f.FaceNodes) {
				bFaces[i].insert(faceNodeIndex);
			};
		};	
		std::cout<<"N="<<nodes.size()<<", Nb="<<bNodes.size()<<'\n';

		//Build kd-tree search index for boundary nodes
		/*
		const int nD = 2;
		int n = 0;
		std::unordered_map<int, int> indexMap;
		flann::Matrix<double> bCoords(new double[bNodes.size() * nD], bNodes.size(), nD);
		for (std::set<int>::iterator it = bNodes.begin(); it != bNodes.end(); it++) {
			//if (Displ.find(*it) == Displ.cend()) continue;
			Node& node = grid.nodes[*it];
			bCoords[n][0] = node.P.x;
			bCoords[n][1] = node.P.y;
			indexMap[n] = node.GlobalIndex;
			n++;
		};

		flann::Index< flann::L2<double> > index(bCoords, flann::KDTreeSingleIndexParams(4));
		index.buildIndex();                                                                                                   
		*/		
	
		/*
		// do a knn search for every inner node, using 128 checks
		int nq = iNodes.size();
		std::vector< std::vector<int> > indices;
		std::vector<std::vector<double> > dists;
		flann::Matrix<double> query(new double[nq * nD], nq, nD);

		n = 0;
		for (std::set<int>::iterator it = iNodes.begin(); it != iNodes.end(); it++) {
			Node& node = grid.nodes[*it];
			query[n][0] = node.P.x;
			query[n][1] = node.P.y;
			n++;
		};
	

		int knn = 50; //number of neighbours used to compute displacement
		index.knnSearch(query, indices, dists, knn, flann::SearchParams(128));
		*/
	
		//Iterate through all boundary nodes compute rotation
		int counter = 0;
		for (int nodeIndex : bNodes) {
			counter++;		
			Node& node = grid.localNodes[nodeIndex];
			Vector Displacement = dR[nodeIndex];
			Vector new_normal;
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
					new_normal += n;
				} else {
					new_normal -= n;
				};			
				normal += f.FaceNormal / f.FaceNormal.mod();
				k++;
			};				
			Vector r = new_normal;
			new_normal /= new_normal.mod();		
			normal /= k;
			normal /= normal.mod();
	//		std::cout<<"Rotation : "<<counter<<" normal" << normal.mod() << "\n";
			RotationMatrix M = CalcRotation(normal, new_normal);
			R[nodeIndex] = M;
			b[nodeIndex] = node.P + Displacement - M * node.P;			
		};

		//Compute displacement for inner nodes
		//n = 0;
		counter = 0;
		double percent = 0;
		for (int nodeIndex : iNodes) {
			counter++;		
			//if (1.0 * counter / iNodes.size() >=  percent) {		
			//	std::cout<<"Interpolation : "<<percent * 100<<"% complete\n";
			//	percent += 0.01;
			//};
			Node& ni = grid.localNodes[nodeIndex];
			std::vector<Vector> displacements;
			dR[nodeIndex] = Vector(0,0,0);
			double sum = 0;		

			//knn-boosted version
			/*
			for (int i=0; i<indices[n].size(); i++) {			
				Node& nb = grid.nodes[indexMap[indices[n][i]]];
				Vector b = *(Vector *)nb.Data["b"];
				RotationMatrix M = *(RotationMatrix *)nb.Data["R"];			
				Vector dr = ni.P - nb.P;
				double w = W(dr);						
				Vector displ = (M * ni.P + b - ni.P);			
				*(Vector *)ni.Data["dR"] += w * displ;			
				sum += w;
			}					
			*(Vector *)ni.Data["dR"] /= sum;
			n++;
			*/
			//Old version
		
			for (int movingNodeIndex : bNodes) {			
				Node& nb = grid.localNodes[movingNodeIndex];
				Vector bb = b[movingNodeIndex];
				RotationMatrix M = R[movingNodeIndex];			
				Vector dr = ni.P - nb.P;
				double w = 1; //W(dr);	TO DO					
				Vector displ = (M * ni.P + bb - ni.P);			
				dR[nodeIndex] += w * displ;			
				sum += w;
			}					
			dR[nodeIndex] /= sum;		
		};

		//Change grid according to displacements	
		//Nodes first
		for (int i = 0;i < grid.localNodes.size(); i++) {
			grid.localNodes[i].P += dR[i];
		};

		//Recalculate Face normals and other grid geometric properties
		grid.UpdateGeometricProperties();
	};
};

#endif