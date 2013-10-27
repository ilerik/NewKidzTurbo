#ifndef TURBO_TURBULENTMODELS_GARBARUKTURBULENTMODEL
#define TURBO_TURBULENTMODELS_GARBARUKTURBULENTMODEL

#include "grid.h"
#include "model.h"

template<class RiemannSolver>
class GarbarukTurbulentPlateModel : public Model<RiemannSolver> {
private:
	DistributedEntityManager<RealValue> TurbulentViscosity;

	//friction velocities and displacement thickness on the wall
	std::map<int, double> DisplacementThickness;
	std::map<int, double> FrictionVelocity;

	//Get local Reinolds number at during cell
	double GetReinolds(Cell& c)
	{
		return U[c.GlobalIndex].ro*U_inf*(c.CellCenter.x - x0)/medium.Viscosity;
	};

	//Compute displacement thickness and friction velocity on wall face
	void ComputeWallVariables1(int wallFaceInd)
	{
		//layer cells must be sorted in ascending order by distance
		std::vector<int> CellSequence = _wallFaces[wallFaceInd].layerCells;
		
		//write required function
		function f(1);
		f.x[0] = 0;		//values on wall
		f.value[0] = 0;

		for(int i=0; i<CellSequence.size(); i++)
		{
			//break condition for monotonic velocity profile
			if(U[CellSequence[i]].rou/U[CellSequence[i]].ro > 0.99*U_inf) break;
			
			//fill in our function
			f.x.push_back(_wallInfo[CellSequence[i]].distance);		
			f.value.push_back(1.0 - U[CellSequence[i]].rou/(U[CellSequence[i]].ro*U_inf));
		};
		f.size = f.x.size();
		//compute integral by trapezium formula
		DisplacementThickness[wallFaceInd] = f.TrapeziumIntegral();
		if(DisplacementThickness[wallFaceInd]==0) throw Exception("Zero displacement Thickness");

		//compute friction velocity
		double tau_wall = medium.Viscosity*(U[CellSequence[0]].rou/U[CellSequence[0]].ro)/_wallInfo[CellSequence[0]].distance;
		FrictionVelocity[wallFaceInd] = sqrt(tau_wall/U[CellSequence[0]].ro);
	};

	//alternative approach
	void ComputeWallVariables(int wallFaceInd)
	{
		//layer cells must be sorted in ascending order by distance in some angular range
		std::vector<int> CellSequence(0);
		for(int i=0; i<_wallFaces[wallFaceInd].layerCells.size(); i++) {
			CellWallInfo& wInfo = _wallInfo[_wallFaces[wallFaceInd].layerCells[i]];
			if(wInfo.angle<0.0000000001) {
				CellSequence.push_back(_wallFaces[wallFaceInd].layerCells[i]);
			};
		};
		
		//write required function
		function f(1);
		f.x[0] = 0;		//values on wall
		f.value[0] = 1;

		//Ќаходим толщину смещени€ из максимальной скорости в срезе!!! ћаксимальан€ скорость может н асовпадать с первым локальным максимумом от стенке
		double Umax = 0;
		int max_cell = 0;
		for(int i=0; i<CellSequence.size(); i++)
		{
			if(U[CellSequence[i]].rou/U[CellSequence[i]].ro > Umax)
			{
				Umax = U[CellSequence[i]].rou/U[CellSequence[i]].ro;
				max_cell++;
			}else break;
		};

		for(int i=0; i<max_cell; i++)
		{
			//fill in our function
			f.x.push_back(_wallInfo[CellSequence[i]].distance);		
			f.value.push_back(1.0 - U[CellSequence[i]].rou/(U[CellSequence[i]].ro*Umax));
		};
		f.size = f.x.size();
		//compute integral by trapezium formula
		DisplacementThickness[wallFaceInd] = f.TrapeziumIntegral();
		if(DisplacementThickness[wallFaceInd]==0) throw Exception("Zero displacement Thickness");

		//compute friction velocity
		double tau_wall = medium.Viscosity*(U[CellSequence[0]].rou/U[CellSequence[0]].ro)/_wallInfo[CellSequence[0]].distance;
		FrictionVelocity[wallFaceInd] = sqrt(tau_wall/U[CellSequence[0]].ro);
	};

	//Compute displacement thickness and friction velocities on all wall
	void ComputeWallVariables()
	{
		//compute displacement thickness and friction velocity in each wall faces and write corresponding maps
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
			ComputeWallVariables(it->first);
	};

	//Compute turbulent viscosity by Garbaruk-Lapin-Strelets model in a cell
	double ComputeTurbViscosityGLS(Cell& c)
	{
		if(GetReinolds(c)<Re_cr) return 0;	//for laminar area
		double y = c.CellCenter.y;  //Wall distance (!!! BRUT IMPLEMENTATION !!!)
		//double y = _wallInfo[c.GlobalIndex].distance;

		//Get displacement thickness and friction velocity at appropriate wall face
		int wallFaceInd = _wallInfo[c.GlobalIndex].wallFaceIndex;
		double Delta = DisplacementThickness[wallFaceInd];
		double U_fr = FrictionVelocity[wallFaceInd];

		//compute all other functions
		double gamma = 1.0/(1 + 5.5*pow((y/Delta), 6));  //Hlebanov function
		double K = 0.41;	//Karman constant
		double A = 12.0;	//for Van Driest function
		double Dvd = pow(1 - exp((-1)*U_fr*y*U[c.GlobalIndex].ro/(medium.Viscosity*A)), 3);	//Van Driest function
		
		return K*U_fr*min(y*Dvd, Delta*gamma);
	};
	
	//and in all cells
	void ComputeTurbViscosityGLS()
	{
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++)
			TurbulentViscosity[cells[i]->GlobalIndex].value = ComputeTurbViscosityGLS(*cells[i]);
	};

	//interpolate turbulent viscosity from cells to face
	double InterpolateTurbulentViscoity(Face& f)
	{
		if(f.FaceCenter.y==0)	return 0;	//for down boarder !!!BRUT IMPLEMENTATION!!!

		//if f is external then turbulent viscosity at face equals to internal one
		double TV_left = TurbulentViscosity[f.FaceCell_1].value;
		if(f.isExternal) return TV_left;

		//for inner faces do interpolation
		double TV_right = TurbulentViscosity[f.FaceCell_2].value;
		//for orthogonal grid
		double r = (_grid.cells[f.FaceCell_1].CellCenter - _grid.cells[f.FaceCell_2].CellCenter).mod();
		double r_l = (f.FaceCenter - _grid.cells[f.FaceCell_1].CellCenter).mod()/r;
		double r_r = (f.FaceCenter - _grid.cells[f.FaceCell_2].CellCenter).mod()/r;

		return r_l*TV_right + r_r*TV_left;
	};

public:
	double U_inf, Re_cr, x0;	//inlet flow velosity, critical Reinolds number and x coordinate of plate border

	//constructor
	GarbarukTurbulentPlateModel(double _U_inf, double _Re_cr, double _x0): U_inf(_U_inf), Re_cr(_Re_cr), x0(_x0) {};

	//bind grid with initializing of turbulent viscocity distribution
	void BindGrid(Grid& grid) {
		//Grid
		_grid = grid;

		//Allocate memory for data over each cell	
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for (int i = 0; i<cells.size(); i++){
			U.add(ConservativeVariables(cells[i]->GlobalIndex));
			TurbulentViscosity.add(RealValue(cells[i]->GlobalIndex));
		};
	};

	//may be not usefull
	double GetTurbulentViscosity(Cell& c)
	{
		return TurbulentViscosity[c.GlobalIndex].value;
	};

	//time step
	void Step() {		
		stepInfo.TimeStep = std::numeric_limits<double>::max();					
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();	
		fluxes.clear();
		for (int i = 0; i<faces.size(); i++) {
			fluxes[faces[i]->GlobalIndex] = std::vector<double>(5,0);
			for (int k = 0; k<5; k++) fluxes[faces[i]->GlobalIndex][k] = 0;
		};
		vfluxes.clear();
		for (int i = 0; i<faces.size(); i++) {
			vfluxes[faces[i]->GlobalIndex] = std::vector<double>(5,0);
			for (int k = 0; k<5; k++) vfluxes[faces[i]->GlobalIndex][k] = 0;
		};

		//Compute gradients of conservative variables for each face
		if (IsGradientRequired) 
			ComputeGradients();	
		
		//compute friction velocities and displacement thickness for all plate
		ComputeWallVariables();

		//compute turbulent viscosity in each cells
		ComputeTurbViscosityGLS();

		//Compute convective fluxes
		//info.TimeStep = 5e-8;
		ComputeConvectiveFluxes();

		//Compute viscous fluxes
		if (IsViscousFluxesRequired)
			ComputeViscousTurbFluxes();

		//Apply CFL condition
		stepInfo.TimeStep *= CFL;
		stepInfo.Residual.resize(5,0);
		for (int i = 0; i<5; i++) stepInfo.Residual[i] = 0;		

		//Distribute fluxes to cells		
		std::vector<double> sumFlux(5,0);		
		for (int i = 0; i<cells.size(); i++) {
			Cell c = *cells[i];
			//std::vector<double> sumFlux2 = sumFluxCell[c.GlobalIndex];

			//For each face			
			for (int k = 0; k<5; k++) sumFlux[k] = 0;

			if (c.Faces.size() != 4) throw Exception("Bu");
			for (int j = 0; j<c.Faces.size(); j++) {
				Face& f = _grid.faces[c.Faces[j]];
				std::vector<double>& flux = fluxes[f.GlobalIndex];					
				std::vector<double>& vflux = vfluxes[f.GlobalIndex];					
				
				//Consider direction
				if (c.GlobalIndex == f.FaceCell_1) {
					sumFlux -= flux + vflux;
				} else {
					sumFlux += flux + vflux;
				};
			};

			//Add sum flux to cell values
			sumFlux = stepInfo.TimeStep * sumFlux / c.CellVolume;
			stepInfo.Residual += abs(sumFlux);

			U[c.GlobalIndex].ro += sumFlux[0];
			U[c.GlobalIndex].rou += sumFlux[1];
			U[c.GlobalIndex].rov += sumFlux[2];
			U[c.GlobalIndex].row += sumFlux[3];
			U[c.GlobalIndex].roE += sumFlux[4];
		};
		
		totalTime += stepInfo.TimeStep;
	};

	//Compute effective viscous flux for each faces
	void ComputeViscousTurbFluxes() {
		//return;
		double max_Du_Dy = 0;

		//Obtain all local faces
		std::vector<Face*> faces = _grid.faces.getLocalNodes();				
		std::vector<double> sumFlux(5,0);

		//turbulent Prandtl number
		double Pr_t = 1.0;

		//Compute viscous flux for each face
		double viscosity = medium.Viscosity;	
		double s_viscosity, effective_viscosity;		//second and effective viscosity
		double turb_viscosity;			//!!!THAT IS TURB VISCOSITY DIVIDED BY DENSITY (KINETIC VISCOSITY)
		double turb_heatcondcoeff;

		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& f = *faces[i];

			if ((f.FaceCenter.y == 0) && (f.FaceCenter.x > 0.96) && (f.FaceCenter.x < 1.01)) {
				std::cout<<f.GlobalIndex<<std::endl;
			};

			std::vector<double> vflux(5, 0);
			for (int k = 0; k<5; k++) vflux[k] = 0;

			std::vector<Vector> grad = gradsFaces[f.GlobalIndex];
			ConservativeVariables UL = U[f.FaceCell_1]; 
            ConservativeVariables UR; 
            if (f.isExternal) { 
                UR = GetDummyRealValues(UL, f); 
            } else { 
                UR = U[f.FaceCell_2]; 
            };
			
			//Here we compute actual flux

			//required variables on face
			double rho = 0.5*(UL[0] + UR[0]);
			double u = 0.5*(UL[1] + UR[1])/rho;
			double v = 0.5*(UL[2] + UR[2])/rho;
			double w = 0.5*(UL[3] + UR[3])/rho;
			double E = 0.5*(UL[4] + UR[4])/rho;

			//compute turbulent viscosity at face f
			turb_viscosity = InterpolateTurbulentViscoity(f);
			effective_viscosity = rho*turb_viscosity + viscosity;
			s_viscosity = -2.0/3.0 * effective_viscosity;

			//TO DO CHECK
			turb_heatcondcoeff = rho*turb_viscosity*medium.Gamma*medium.Cv/Pr_t;
				
			//gradients of required variables
			Vector du, dv, dw, dE, dT;
			du.x = (grad[1].x - u*grad[0].x)/rho;
			du.y = (grad[1].y - u*grad[0].y)/rho;
			du.z = (grad[1].z - u*grad[0].z)/rho;

			dv.x = (grad[2].x - v*grad[0].x)/rho;
			dv.y = (grad[2].y - v*grad[0].y)/rho;
			dv.z = (grad[2].z - v*grad[0].z)/rho;

			dw.x = (grad[3].x - w*grad[0].x)/rho;
			dw.y = (grad[3].y - w*grad[0].y)/rho;
			dw.z = (grad[3].z - w*grad[0].z)/rho;

			dE.x = (grad[4].x - E*grad[0].x)/rho;
			dE.y = (grad[4].y - E*grad[0].y)/rho;
			dE.z = (grad[4].z - E*grad[0].z)/rho;

			double R = (medium.Gamma - 1.0) * medium.Cv;		//specific gas constant
			double K = medium.ThermalConductivity + turb_heatcondcoeff;		//effective heat conduction coefficient
			dT.x = (medium.Gamma - 1)*(dE.x - u*du.x - v*dv.x - w*dw.x)/R;
			dT.y = (medium.Gamma - 1)*(dE.y - u*du.y - v*dv.y - w*dw.y)/R;
			dT.z = (medium.Gamma - 1)*(dE.z - u*du.z - v*dv.z - w*dw.z)/R;

			//
			if (abs(du.y) > max_Du_Dy) max_Du_Dy = abs(du.y);

			//stress elements
			double tau_diagonal = s_viscosity*(du.x + dv.y + dw.z);
			double tau_xx = tau_diagonal + 2*effective_viscosity*du.x;
			double tau_yy = tau_diagonal + 2*effective_viscosity*dv.y;
			double tau_zz = tau_diagonal + 2*effective_viscosity*dw.z;
			double tau_xy = effective_viscosity*(du.y + dv.x);
			double tau_xz = effective_viscosity*(du.z + dw.x);
			double tau_yz = effective_viscosity*(dv.z + dw.y);

			//work of viscous stresses and heat conduction
			Vector Thetta;
			Thetta.x = u*tau_xx + v*tau_xy + w*tau_xz + K*dT.x;
			Thetta.y = u*tau_xy + v*tau_yy + w*tau_yz + K*dT.y;
			Thetta.z = u*tau_xz + v*tau_yz + w*tau_zz + K*dT.z;

			//TO DO Remove fix for symmetry boundary
			if (f.BCMarker == 3) {
				tau_xy = 0;
				tau_xz = 0;
				Thetta.y = 0;
			};

			//viscous fluxes
 			vflux[1] = f.FaceNormal.x*tau_xx + f.FaceNormal.y*tau_xy + f.FaceNormal.z*tau_xz;			
			vflux[2] = f.FaceNormal.x*tau_xy + f.FaceNormal.y*tau_yy + f.FaceNormal.z*tau_yz;			
			vflux[3] = f.FaceNormal.x*tau_xz + f.FaceNormal.y*tau_yz + f.FaceNormal.z*tau_zz;			
			vflux[4] = f.FaceNormal*Thetta;
			
			//
			for (int k = 0; k<5; k++) vflux[k] *= f.FaceSquare;

			//Add viscous fluxes		
			sumFlux = abs(vflux);
			std::vector<double> cflux = fluxes[f.GlobalIndex];
			//fluxes[f.GlobalIndex] -= vflux; //“ут по идее он должен быть умноже на площадь грани так что обрати внимание на свои формулы

			//TO DO remove
			vfluxes[f.GlobalIndex] = vflux;
		};

		std::cout<<"VFlux[1] = " << sumFlux[1] << std::endl;
		std::cout<<"Max du/dy = " << max_Du_Dy << std::endl;
	};

	//add new variable - turbulent viscocity
	void SaveToTechPlot(std::string fname)
 {
		std::ofstream ofs(fname);
		//TO DO unify		
		//Header
		ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\", \"distance\", \"t_viscosity\"\n";
		ofs<<"ZONE T=\"D\"\n";
		ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ET=QUADRILATERAL\n";
		ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED)\n";

		//Map all node global indexes to natural numbers
		std::map<int,int> toNaturalIndex;
		std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
		int counter = 1;
		for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;


		//Access local nodes, faces, cells and flow data
		std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
		std::vector<Face*> faces = _grid.faces.getLocalNodes();
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		//Nodes coordinates
		//X
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.x<<"\n";
		};

		//Y
		for (int i = 0; i<nodes.size(); i++) {
			ofs<<nodes[i]->P.y<<"\n";
		};

		////Solution data
		////Rho
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].ro<<"\n";
		};
				
		//u
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].rou/U[idx].ro<<"\n";
		};
		//v
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].rov/U[idx].ro<<"\n";		
		};
		//w
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].row/U[idx].ro<<"\n";
		};
		//E
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<U[idx].roE/U[idx].ro<<"\n";
		};
		//T
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double T = (E - (vx*vx+vy*vy+vz*vz)/2.0) / (medium.Cv);
			ofs<<T<<"\n";
		};
		//P
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double P = (medium.Gamma-1.0) * ro * (E - (vx*vx+vy*vy+vz*vz)/2.0);
			ofs<<P<<"\n";
		};		

		//distance
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			if (_wallInfo.size() != 0) {
				ofs<<_wallInfo[idx].distance<<"\n";
			} else {
				ofs<<0<<"\n";
			};
		};

		//turbulent viscosity
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<TurbulentViscosity[idx].value<<"\n";
		};

		//Connectivity list for each cell
		for (int i = 0; i<cells.size(); i++) {
			ofs<<toNaturalIndex[cells[i]->Nodes[0]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[1]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[2]]<<" "
				<<toNaturalIndex[cells[i]->Nodes[3]]<<"\n";
		};

		ofs.close();
		return;

	};

	//function for testing turb viscosity calculation
	void Test()
	{
		//ComputeWallDistances();
		ComputeWallVariables();
		ComputeTurbViscosityGLS();
		SaveToTechPlot("turb_viscosity.dat");
	};
};

#endif