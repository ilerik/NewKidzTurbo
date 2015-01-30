#ifndef TURBO_ALEmethod_ALEMethod
#define TURBO_ALEmethod_ALEMethod

#include "grid.h"
#include "utilityfunctions.h"
#include <unordered_set>

//Define ALE integration method
class ALEMethod {	
	Grid* _grid;

	// map faceIndex to [0,1] indicator value. 1.0 - F
	double (*_indicatorFunction)(int); 
public:
	//Set pointer to indicator function

	//Set reference to grid
	void SetGrid(Grid& grid) {
		_grid = &grid;
	};

	//Define availible types of ALE motion
	//Method type
	enum class ALEMotionType {
		PureEulerian,
		PureLagrangian,
		ALEMaterialInterfaces
	} ALEMotionType;

	//Velocity of faces and nodes for ALE step
	std::map<int, Vector> nodesVelocity;
	std::map<int, Vector> facesVelocity;
	std::map<int, double> facesPressure; 

	//Lists of nodes with defined motion and free nodes
	std::unordered_set<int> movingNodes;
	std::unordered_set<int> freeNodes;

	//List of adjacent faces for each node
	std::map<int, std::vector<int>> adjacentFaces;

	//Compute swept volume
	Vector ComputeFaceVelocityByNodes(Face& f) {
		Vector velocity = Vector(0,0,0);
		if (f.CGNSType == NODE) {
			velocity = nodesVelocity[f.FaceNodes[0]];
			return velocity;
		};
		if (f.CGNSType == BAR_2) {
			Vector v1 = nodesVelocity[f.FaceNodes[0]];
			Vector v2 = nodesVelocity[f.FaceNodes[1]];
			velocity = 0.5 * (v1+v2);
			return velocity;
		};
		throw new Exception("Unsupported element type");
		return velocity;
	};

	//ALE residual correction procedure
	void Remap(int nVariables, std::vector<double>& residual, std::vector<double>& cellValues, double dt) {
		//Compute grid motion fluxes through each face
		for (int cellInd = 0; cellInd < _grid->nCellsLocal; cellInd++) {
			Cell* c = _grid->localCells[cellInd];
			for (int faceInd : c->Faces) {
				Face& f = _grid->localFaces[faceInd];					

				//Compute volume swept by face TO DO generalize
				double sweptVolume = 0;
				Vector avgVelocity = Vector(0,0,0);
				for (int nodeInd : f.FaceNodes) {
					avgVelocity += nodesVelocity[nodeInd];
				};
				avgVelocity /= f.FaceNodes.size();
				sweptVolume = avgVelocity.mod() * dt;

				//Distribute moving face flux				
				for (int i = 0; i<nVariables; i++) { 
					residual += cellValues * sweptVolume;
				};

				//Vector uR = ;
				//double roCell = Values[cellInd * nVariables + 0];					
				//double A = -1.0 * f.FaceSquare;
				//if (f.FaceCell_1 != c->GlobalIndex) {
				//	A *= -1; //Reverse flow if normal directed inwards
				//	//roCell = GetCellValues(f.FaceCell_2)[0];
				//};
				////A *= FaceFluxes[f.GlobalIndex][0] / (roCell * un);

				//for (int j = 0; j<nVariables; j++) {
				//	double dR =  -Values[cellInd * nVariables + j] * un * A;
				//	residual[cellInd * nVariables + j] += dR; // / _grid.localCells[cellInd]->CellVolume;
				//};
			};				
		};
	};

	//Compute face velocity
	Vector ComputeFaceVelocity(GasModel* gasModelLeft, const GasModel::ConservativeVariables& UL, GasModel* gasModelRight, const GasModel::ConservativeVariables& UR, const Face& f, double ALEindicator) {
		Vector faceVelocity;
		if (ALEindicator == 0) {
			faceVelocity = Vector(0,0,0);
			return faceVelocity;
		} else {
			Vector faceNormalVelocity;
			Vector faceTangentialVelocity;

			//The Harten, Lax, and van Leer with contact restoration (HLLC) Riemann solver
			//Left and right states
			double roL = UL.ro;
			Vector vL = Vector(UL.rou / roL, UL.rov / roL, UL.row / roL);
			double uL = vL * f.FaceNormal;
			double pL = 0;
			double cL = 0;
			double GrL = 0;	
			assert(roL > 0);
			gasModelLeft->GetPressureAndSoundSpeed(UL, pL, cL, GrL);				
			double phiL = cL*cL - GrL * pL / roL;

			double roR = UR.ro;
			Vector vR = Vector(UR.rou / roR, UR.rov / roR, UR.row / roR);
			double uR = vR * f.FaceNormal;
			double pR = 0;
			double cR = 0;
			double GrR = 0;
			double GammaR = 0;				
			assert(roR > 0);
			gasModelRight->GetPressureAndSoundSpeed(UR, pR, cR, GrR);			
			double phiR = cR*cR - GrR * pR / roR;

			//Generalized Roe averages (according to Hu et al)
			double roLRoot = sqrt(roL);
			double roRRoot = sqrt(roR);
			double roRoe = roLRoot * roRRoot; //Roe averaged density
			double uRoe = takeRoeAverage(roLRoot, roRRoot, uL, uR); //Roe averaged velocity (Hu et al (17))
			double pOverRoRoe = takeRoeAverage(roLRoot, roRRoot, pL/roL, pR/roR) + 0.5*pow((uR - uL) / (roRRoot + roLRoot), 2.0); //(Hu et al (18))		
			double phiRoe = takeRoeAverage(roLRoot, roRRoot, phiL, phiR); //Roe averaged phi (Hu et al (25))
			double GrRoe = takeRoeAverage(roLRoot, roRRoot, GrL, GrR); //Roe averaged Gruneisen coefficient (Hu et al (25))
			double cRoe = phiRoe + GrRoe * pOverRoRoe; //Roe averaged sound speed (Hu et al (16))
			cRoe = sqrt(cRoe);

			//Wave speeds for two waves (Hu et al (12))
			double SL = min(uL - cL, uRoe - cRoe); //bl
			double SR = max(uR + cR, uRoe + cRoe); //br

			//The HLLC Flux for Euler equations (Toro 10.4.2, p324)
			double SStar = (pR - pL) + roL*uL*(SL - uL) - roR*uR*(SR - uR); //Toro 10.37
			SStar /=  roL*(SL - uL) - roR*(SR - uR); // Speed of middle wave

			//Now obtain face velocity
			faceNormalVelocity = f.FaceNormal * SStar;
			Vector tangVelocityL = vL - uL * f.FaceNormal;
			Vector tangVelocityR = vR - uR * f.FaceNormal;
			faceTangentialVelocity = (tangVelocityL + tangVelocityR) / 2.0;
			faceVelocity = ALEindicator * (faceNormalVelocity + faceTangentialVelocity);

			double PStarL = pL + roL * (SL - uL)*(SStar - uL); //Pressure estimate for left star region Toro 10.36
			double PStarR = pR + roR * (SR - uR)*(SStar - uR); //Pressure estimate for right star region Toro 10.36

			return faceVelocity;
		};
	};

	//Compute velocities
	void ComputeMovingNodesVelocities() {
		movingNodes.clear();

		//For each node refresh list of adjacent faces with computed velocities
		adjacentFaces.clear();		
		for (Face& f : _grid->localFaces) {			
			//If one of the faces is moving then it's moving node
			bool isMovingNode = false;
			if (facesVelocity.find(f.GlobalIndex) != facesVelocity.end()) {
				isMovingNode = true;
			};

			//If face is boundary then it's nodes also considered moving
			if (_grid->IsBoundaryFace(f)) {
				isMovingNode = true;
			};
			
			//Store adjacent faces for each node
			for (int node : f.FaceNodes) {				
				adjacentFaces[node].push_back(f.GlobalIndex);
				if (isMovingNode) movingNodes.insert(node);
			};
		};

		//Determine free nodes
		//TO DO make optimal
		freeNodes.clear();
		for (Node& n : _grid->localNodes) {
			if (movingNodes.find(n.GlobalIndex) == std::end(movingNodes)) freeNodes.insert(n.GlobalIndex);
		};

		//For each moving node compute velocity (node is moving if any adjacent face is moving)
		for (int nodeIndex : movingNodes) {	
			Node& n = _grid->localNodes[nodeIndex];
			Vector nodeVelocity; //Resulting node velocity
			double sumFacesSquare = 0;

			//Matrix consisting of face normals
			std::vector<double> faceNormals;
			//Vector of face normal velocities
			std::vector<double> faceVelocities;

			std::vector<double> velocities;
			std::vector<Vector> normals;
			std::vector<double> weights;
			
			//LLS approximation
			for (int faceIndex : adjacentFaces[n.GlobalIndex]) {
				Face& face = _grid->localFaces[faceIndex];
				Vector faceVelocity = facesVelocity[faceIndex];
				velocities.push_back(faceVelocity * face.FaceNormal);
				normals.push_back(face.FaceNormal);
				weights.push_back(face.FaceSquare);				
			};
			nodeVelocity = ComputeVelocityByPoints(_grid->gridInfo.CellDimensions, velocities, normals, weights);

			//Store node velocity
			nodesVelocity[n.GlobalIndex] = nodeVelocity;
		};	
	};

	//Compute free nodes velocities
	void ComputeFreeNodesVelocities() {
		//TO DO 1D implementation temporary
		std::vector<double> coordinates;
		std::map<double, Vector> velocitiesByCoordinate;

		for (int nodeIndex : movingNodes) {
			double coord = _grid->localNodes[nodeIndex].P.x;
			coordinates.push_back(coord);
			velocitiesByCoordinate[coord] = nodesVelocity[nodeIndex];
		};
		
		std::sort(coordinates.begin(), coordinates.end());

		for (int nodeIndex : freeNodes) {
			Vector nodeVelocity; //Resulting node velocity
			double x = _grid->localNodes[nodeIndex].P.x;
			for (int i = 1; i < coordinates.size(); i++) {
				double xL = coordinates[i-1];
				double xR = coordinates[i];
				Vector vL = velocitiesByCoordinate[xL];
				Vector vR = velocitiesByCoordinate[xR];
				if ((xR >= x) && (xL <= x)) {
					nodeVelocity = vL + (x - xL) * (vR - vL) / (xR - xL);
				};
			};
			nodesVelocity[nodeIndex] = nodeVelocity;
		};
	};

	//Move mesh
	void MoveMesh(double timestep) {
		for (std::pair<int, Vector> pair : nodesVelocity) {
			int nodeInd = pair.first;
			Vector v = pair.second;
			_grid->localNodes[nodeInd].P += v * timestep;
		};
	};

	//Remeshing procedure
	void Remesh() {

	};

};

#endif