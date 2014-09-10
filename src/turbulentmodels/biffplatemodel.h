#ifndef TURBO_TURBULENTMODELS_BIFFURCATIONTURBULENTMODEL
#define TURBO_TURBULENTMODELS_BIFFURCATIONTURBULENTMODEL

#include "grid.h"
#include "model.h"

//functions and parameters for Levenberg–Marquardt solver
function BM_FUN_vel_der;
function BM_FUN_t_visc;
function BM_FUN_t_dif;
function BM_FUN_t_dif_der;
function BM_FUN_lh_vel_der;
double BM_PAR_eps, BM_PAR_alpha, BM_PAR_c;

//compute residual and write it in x array
void ComputeWriteRes(const alglib::real_1d_array& x, alglib::real_1d_array &f, void *ptr)
{
	std::vector<double> grid = BM_FUN_vel_der.x;

	//left bounday condition
	f[0] = x[0];
	for(int i=1; i<f.length() - 1; i++)
	{
		//left and right space step
		double hl = grid[i] - grid[i-1];
		double hr = grid[i+1] - grid[i];
		double h = hr + hl;

		//compute equation coefficients
		double am = 2*(BM_PAR_eps + BM_FUN_t_dif.value[i])*BM_FUN_lh_vel_der.value[i] - BM_FUN_t_dif_der.value[i]*BM_FUN_lh_vel_der.value[i]*hr;
		am *= BM_PAR_eps*BM_FUN_t_visc.value[i]/(hl*h);

		double ap = 2*(BM_PAR_eps + BM_FUN_t_dif.value[i])*BM_FUN_lh_vel_der.value[i] + BM_FUN_t_dif_der.value[i]*BM_FUN_lh_vel_der.value[i]*hl;
		ap *= BM_PAR_eps*BM_FUN_t_visc.value[i]/(hr*h);

		double a = BM_FUN_t_dif_der.value[i]*(hr - hl)/(hr*hl) - 2*(BM_PAR_eps + BM_FUN_t_dif.value[i])/(hr*hl) - BM_PAR_eps*BM_PAR_alpha;
		a *= BM_PAR_eps*BM_FUN_lh_vel_der.value[i];
		a += BM_PAR_c*BM_PAR_c*BM_FUN_vel_der.value[i];
		a *= BM_FUN_t_visc.value[i];

		//compute a components of optimized function
		f[i] = am*x[i-1] + a*x[i] + ap*x[i+1] + BM_PAR_c*BM_PAR_c*x[i]*x[i];
	};
	//right boundary condition
	f[f.length() - 1] = x[f.length() - 1];
	return; 
};
//alternative target function (the same one with longer computation time)
void ComputeWriteRes2(const alglib::real_1d_array& x, alglib::real_1d_array &f, void *ptr)
{
	std::vector<double> grid = BM_FUN_vel_der.x;

	//derivatives of x
	function x_der(x.length());
	x_der.x = grid;
	for(int i=0; i<x_der.size; i++) x_der.value[i] = x[i];
	function x_sec_der = x_der.DerivateSec();
	x_der = x_der.Derivate();

	//left boundary condition
	f[0] = x[0];
	for(int i=1; i<f.length() - 1; i++)
	{
		double f1 = BM_FUN_t_dif_der.value[i]*x_der.value[i] + BM_FUN_t_dif.value[i]*x_sec_der.value[i];
		f1 *= (-1.0)*BM_PAR_eps*BM_FUN_lh_vel_der.value[i]*BM_FUN_t_visc.value[i];

		double f2 = BM_FUN_vel_der.value[i]*BM_FUN_t_visc.value[i] + x[i];
		f2 *= (-1.0)*BM_PAR_c*BM_PAR_c*x[i];

		double f3 = BM_FUN_t_visc.value[i]*BM_FUN_lh_vel_der.value[i]*BM_PAR_eps*BM_PAR_eps;
		f3 *= BM_PAR_alpha*x[i];
		
		//compute a components of optimized function
		f[i] = f1 + f2 + f3;
	};
	//right boundary condition
	f[f.length() - 1] = x[f.length() - 1];
	return; 
};


template<class RiemannSolver>
class BiffTurbulentPlateModel : public Model<RiemannSolver> {
private:
	//Reinolds correlations
	std::map<int, double> xy_turb_stressesC;	//for cells
	std::map<int, double> xy_turb_stressesF;	//for faces

	//Dinamic Re numbers
	std::map<int, double> Re_din;

	//main characteristics
	double U_inf;
	double h;

	//Get local Reinolds number at during cell
	double GetReinolds(Cell& c)
	{
		return U[c.GlobalIndex].ro*U_inf*(c.CellCenter.x - x0)/medium.Viscosity;
	};

	//turbulent viscosity and diffusion coefficients
	function ComputeTurbViscosity(const std::vector<int>& CellSequence, const function& velosity)
	{
		if(velosity.size==0) throw Exception("velocity function has zero size\n");

		//turbulent viscosity function
		function res(velosity.size);
		res.value[0] = 0;
		res.x[0] = 0;

		//wall face for this cells layer
		double Delta = ComputeDisplacementThickness(velosity);
		double tau_wall = medium.Viscosity*velosity.value[1]/velosity.x[1];		//first order
		double U_fr = sqrt(tau_wall/U[CellSequence[0]].ro);
		
		for(int i=1; i<res.size; i++)
		{
			//compute all other functions
			double gamma = 1.0/(1 + 5.5*pow((velosity.x[i]/h), 6));  //Hlebanov function
			double K = 0.41;	//Karman constant
			double A = 12.0;	//for Van Driest function
			double Dvd = pow(1 - exp((-1)*U_fr*velosity.x[i]*U[CellSequence[i-1]].ro/(medium.Viscosity*A)), 3);	//Van Driest function
		
			res.value[i] = K*U_fr*min(velosity.x[i]*Dvd, Delta*gamma);
			res.x[i] = velosity.x[i];
		};
		return res;
	};
	function ComputeTurbDiffusion(const std::vector<int>& CellSequence, const function& velocity)
	{
		return ComputeTurbViscosity(CellSequence, velocity);
	};
	function ComputeLeeHarshaVelocity(const function& velocity)
	{
		return velocity;
	}

	//LEVENBERG-MARQUARDT METHOD for solving Reinolds stress equation
	function ComputeReinoldsStresses()
	{
		//initial approach
		alglib::real_1d_array tau_un;	//that's not stresses but pulses product
		tau_un.setlength(BM_FUN_vel_der.size);
		tau_un[0] = 0;
		for(int i=1; i<tau_un.length() - 1; i++)
		{
			tau_un[i] = (-1.0)*BM_FUN_t_visc.value[i]*BM_FUN_vel_der.value[i];
			//tau_un[i] = -1.0; //alternative initial approach
		};
		tau_un[tau_un.length() - 1] = 0;
		
		//double epsg = 0.0000000001;
		double epsg = 0.000000000000000000001;
		double epsf = 0;
		double epsx = 0;
		alglib::ae_int_t maxits = 0;
		alglib::ae_int_t size = tau_un.length();
		alglib::minlmstate state;
		alglib::minlmreport rep;
		double diffstep = 0.0001;

		alglib::minlmcreatev(size, tau_un, diffstep, state);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
		alglib::minlmoptimize(state, ComputeWriteRes);
		alglib::minlmresults(state, tau_un, rep);

		//printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
		//printf("%d\n", rep.iterationscount);
		//printf("%d\n", rep.nfunc);
		//printf("%s\n", tau_un.tostring(size).c_str()); // EXPECTED: [-3,+3]

		function res = BM_FUN_vel_der;
		for(int i=0; i<res.size; i++)
			res.value[i] = tau_un[i];
		
		return res;
	};

	//create new special grid
	std::vector<int> ComputeGrid(int N)
	{

	};

	//transform dimension variables to dimensionless ones and vice versa
	void TransformDimensionVar()
	{
		BM_FUN_vel_der.value = BM_FUN_vel_der.value*(h/U_inf);
		BM_FUN_vel_der.x = BM_FUN_vel_der.x/h;

		BM_FUN_t_visc.value = BM_FUN_t_visc.value/(U_inf*h);
		BM_FUN_t_visc.x = BM_FUN_t_visc.x/h;

		BM_FUN_t_dif.value = BM_FUN_t_dif.value/(U_inf*h);
		BM_FUN_t_dif.x = BM_FUN_t_dif.x/h;

		BM_FUN_t_dif_der.value = BM_FUN_t_dif_der.value/(U_inf);
		BM_FUN_t_dif_der.x = BM_FUN_t_dif_der.x/h;

		BM_FUN_lh_vel_der.value = BM_FUN_lh_vel_der.value*(h/U_inf);
		BM_FUN_lh_vel_der.x = BM_FUN_lh_vel_der.x/h;

		return;
	};
	void TransformDimensionlessVar()
	{
		BM_FUN_vel_der.value = BM_FUN_vel_der.value*(U_inf/h);
		BM_FUN_vel_der.x = BM_FUN_vel_der.x*h;

		BM_FUN_t_visc.value = BM_FUN_t_visc.value*(U_inf*h);
		BM_FUN_t_visc.x = BM_FUN_t_visc.x*h;

		BM_FUN_t_dif.value = BM_FUN_t_dif.value*(U_inf*h);
		BM_FUN_t_dif.x = BM_FUN_t_dif.x*h;

		BM_FUN_t_dif_der.value = BM_FUN_t_dif_der.value*(U_inf);
		BM_FUN_t_dif_der.x = BM_FUN_t_dif_der.x*h;

		BM_FUN_lh_vel_der.value = BM_FUN_lh_vel_der.value*(U_inf/h);
		BM_FUN_lh_vel_der.x = BM_FUN_lh_vel_der.x*h;

		return;
	};

	//Compute turbulent stresses in the layer of cells over wall face f
	void ComputeTurbStresses(int wallFaceInd)
	{
		//layer cells must be sorted in ascending order by distance in some angular range
		std::vector<int> CellSequence(0);
		U_inf = 0;
		Vector norm = _grid.faces[wallFaceInd].FaceNormal;

		//to compute average density
		function rho;
		//velocity function
		function v(1);
		v.x[0] = 0;
		v.value[0] = 0;
		//layer of cells
		std::vector<int>& layerCells = _wallFaces[wallFaceInd].layerCells;

		//fill CellSequence and velocity function
		for(int i=0; i<layerCells.size(); i++)
		{
			CellWallInfo& wInfo = _wallInfo[layerCells[i]];
			if(wInfo.angle<SoartingAngle) {
				int c_Ind = layerCells[i];
				//tangential velocity in a cell
				double U_tan = GetVelocityX(U[c_Ind]);	//only for Plate!
				if(U_tan > U_inf)	//condition for height of the boundary layer
				{
					CellSequence.push_back(c_Ind);
					v.value.push_back(U_tan);
					v.x.push_back(wInfo.distance*cos(wInfo.angle));
					rho.value.push_back(U[c_Ind].ro);
					rho.x.push_back(wInfo.distance);
					U_inf = U_tan;
				}else break;
			};
		};
		v.size = v.x.size();
		rho.size = rho.x.size();

		//turbulent layer height
		h = v.x[v.size-1];
		double rho_ave;
		if(rho.size>1) rho_ave = rho.TrapeziumIntegral()/(rho.x[rho.size-1] - rho.x[0]);
		else rho_ave = rho.value[0];
		double kinetic_visc = medium.Viscosity/rho_ave;
		_boundaryLayerHeight[wallFaceInd] = h;
		BM_PAR_eps = kinetic_visc/(h*U_inf);
		Re_din[wallFaceInd] = 1.0/BM_PAR_eps;

		//laminar flow condition
		if(BM_PAR_eps>BM_PAR_c/sqrt(BM_PAR_alpha)) return;
		//if(v.size<3) throw Exception("Turbulent boundary layer has not enough grid nodes\n");
		if(v.size<3) return;	//WID çàãëóøêà í àñëó÷àé âûëåòà â ñòðî÷êå âûøå

		//compute usefull functions on general grid
		BM_FUN_t_visc = ComputeTurbViscosity(CellSequence, v);
		//if(wallFaceInd==114) BM_FUN_t_visc.WriteData("t_visc_srez_BM.dat");
		BM_FUN_t_dif = BM_FUN_t_visc;	// !!! only if t_dif = t_visc !!!
		BM_FUN_t_dif_der = BM_FUN_t_dif.Derivate();
		BM_FUN_vel_der = v.Derivate();
		BM_FUN_lh_vel_der = BM_FUN_vel_der;	// !!! only if LeeHersha Velocity = Velocity !!!
		//BM_FUN_t_visc.WriteData("T_visc.dat"); //TO DO DELETE
		
		//transform to dimensionless variables
		TransformDimensionVar();
		function tau_un = ComputeReinoldsStresses();
		//tau_un.WriteData("tau_un.dat");		//TO DO DELETE

		//transform stresses to phisical form
		for(int i = 0; i<CellSequence.size(); i++)
			//xy_turb_stressesC[CellSequence[i]] = U_inf*U_inf*tau_un.value[i+1]*rho.value[i];
			xy_turb_stressesC[CellSequence[i]] = (-1.0)*U_inf*U_inf*tau_un.value[i+1]*rho.value[i];

		return;
	};
	//and in all cells
	void ComputeTurbStresses()
	{
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		//zero initial value in all cells
		for(int i=0; i<cells.size(); i++)
			xy_turb_stressesC[cells[i]->GlobalIndex] = 0;

		//compute stresses over all wall faces
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
			ComputeTurbStresses(it->first);

		//interpolate stresses to all faces
		std::vector<Face*> faces = _grid.faces.getLocalNodes();
		#pragma omp for
		for (int i = 0; i<faces.size(); i++) {
			Face& f = *faces[i];
			double Str_l = xy_turb_stressesC[f.FaceCell_1];
			double r_l = (f.FaceCenter - _grid.cells[f.FaceCell_1].CellCenter).mod();	//WID
			double Str_r, r_r;
			if (f.isExternal) {
				if(IsWallFace(f)) Str_r = (-1.0)*Str_l;
				else Str_r = Str_l;
				r_r = r_l;
			} else {
				Str_r = xy_turb_stressesC[f.FaceCell_2];
				r_r = (f.FaceCenter - _grid.cells[f.FaceCell_2].CellCenter).mod();
			};
			double r = r_l + r_r;	r_l /= r;	r_r /= r;
			xy_turb_stressesF[f.GlobalIndex] = r_l*Str_r + r_r*Str_l;
			xy_turb_stressesF[f.GlobalIndex] = xy_turb_stressesF[f.GlobalIndex];
		};
	};

public:
	double Re_cr;	//critical Reinolds number

	//Compute Reinolds stresses
	void AditionalStepComputations()
	{
		ComputeTurbStresses();
	};

	//get stresses on the face
	double GetTurbViscosity(int FaceInd)
	{
		return xy_turb_stressesF[FaceInd];
	};

	//compute turbulent heat conductivity coefficient
	double GetTurbHeatConductivity(int FaceInd)
	{
		return 0;
		//TO DO use bool parameter of turbulent or laminar model
		double turbulentPrandtl = 0.677;
		double KTurb = medium.Cv * medium.Gamma * GetTurbViscosity(FaceInd) / turbulentPrandtl;
		return KTurb;
	};
	Matrix GetAdditionalStresses(int FaceInd)
	{
		Matrix res(3,3);	//zero matrix
		res[0][1] = xy_turb_stressesF[FaceInd];
		res[1][0] = res[0][1];	//for symmetry
		return res;
	};

	//constructor
	BiffTurbulentPlateModel(double _Re_cr): Re_cr(_Re_cr){};

	BiffTurbulentPlateModel(double _Re_cr, double alpha, double c): Re_cr(_Re_cr)
	{
		BM_PAR_alpha = alpha;
		BM_PAR_c = c;
	};

	//model settings 
	void SetAlpha(double alpha)
	{
		BM_PAR_alpha = alpha;
	};
	void SetC(double c)
	{
		BM_PAR_c = c;
	};

	void SaveToTechPlot(std::string fname)
 {
		std::ofstream ofs(fname);
		//TO DO unify		
		//Header
		ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\", \"distance\", \"Reinolds_stresses\"\n";
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

		//Reinolds stresses
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<xy_turb_stressesC[idx]<<"\n";
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
	void SaveReDinToTechPlot(std::string fname = "Re_din.dat")
	{
		std::ofstream ofs(fname);
		for (std::map<int, FaceWallInfo>::iterator it = _wallFaces.begin(); it != _wallFaces.end(); it++)
		{
			int FaceInd = it->first;
			ofs << _grid.faces[FaceInd].FaceCenter.x << ' ' << Re_din[FaceInd] << '\n';
		};
		ofs.close();
	};
	//??? what is the function?
	void SaveToTecplot25D(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;

		//Access local nodes, faces, cells and flow data
		std::vector<Node*> nodesOriginal = _grid.nodes.getLocalNodes();
		std::vector<Face*> facesOriginal = _grid.faces.getLocalNodes();
		std::vector<Cell*> cellsOriginal = _grid.cells.getLocalNodes();
		std::vector<ConservativeVariables*> data = U.getLocalNodes();

		//Collect everything from Z=0 plane
		std::map<int,int> toNaturalIndex;			
		int counter = 1;
		std::vector<Node*> nodes;
		std::vector<Cell*> cells;
		for (int i = 0; i<nodesOriginal.size(); i++) {
			if (nodesOriginal[i]->P.z < 1e-2) {
				nodes.push_back(nodesOriginal[i]);
				//Map all node global indexes to natural numbers
				toNaturalIndex[nodesOriginal[i]->GlobalIndex] = counter++;
			};
		};

		std::map<int, std::vector<int>> cellPlaneFaceNodes;
		for (int i = 0; i<cellsOriginal.size(); i++) {
			cells.push_back(cellsOriginal[i]);
			//Find face in Z=0 plane
			bool isFound = false;
			Cell& c = *cellsOriginal[i];
			for (int j = 0; j<c.Faces.size(); j++) {
				if (_grid.faces[c.Faces[j]].FaceCenter.z < 1e-2) {
					cellPlaneFaceNodes[c.GlobalIndex] = _grid.faces[c.Faces[j]].FaceNodes;
					isFound = true;
					break;
				};				
			};
			if (!isFound) {
				std::cout<<"Face in plane Z=0 not found\n";
			};
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][0]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][1]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][2]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
			//if (toNaturalIndex.find(cellPlaneFaceNodes[cellsOriginal[i]->GlobalIndex][3]) == toNaturalIndex.end()) std::cout<<"Node not from Z=0 plane\n";
		};

		//TO DO unify		
		//Header
		ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"w\", \"E\", \"T\", \"P\"";
		ofs<<", \"PStagnation\"";
		ofs<<", \"distance\"";
		ofs<<", \"Reinolds_stresses\"";
		ofs<<"\n";
			
		ofs<<"ZONE T=\"D\"\n";
		ofs<<"N=" << nodes.size() << ", E=" << cells.size() <<", F=FEBLOCK, ET=QUADRILATERAL\n";

		ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<", CELLCENTERED";
		ofs<<")\n";

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

		//PStagnation
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double ro = U[idx].ro;
			double vx = U[idx].rou/U[idx].ro;
			double vy = U[idx].rov/U[idx].ro;
			double vz = U[idx].row/U[idx].ro;
			double E = U[idx].roE/U[idx].ro;
			double PStagnation = (medium.Gamma-1.0) * U[idx].roE;
			ofs<<PStagnation<<"\n";
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

		//Reinolds stresses
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<xy_turb_stressesC[idx]<<"\n";
		};

		//Connectivity list for each cell
		for (int i = 0; i<cells.size(); i++) {
			ofs<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][0]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][1]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][2]]<<" "
				<<toNaturalIndex[cellPlaneFaceNodes[cells[i]->GlobalIndex][3]]<<"\n";
		};		

		ofs.close();
		return;
	};

	//function to TESTIFY model
	//write main functions which used to solve Reinolds stress equation
	void WriteArgumentFunctions(int wallFaceInd)
	{
		std::map<int, double> U_top;
		U_top[wallFaceInd] = 0;
		//layer cells must be sorted in ascending order by distance in some angular range
		std::vector<int> CellSequence(0);

		//velocity function
		function v(1);
		v.x[0] = 0;
		v.value[0] = 0;

		//fill in layerCells and velocity function
		for(int i=0; i<_wallFaces[wallFaceInd].layerCells.size(); i++)
		{
			CellWallInfo& wInfo = _wallInfo[_wallFaces[wallFaceInd].layerCells[i]];
			if(wInfo.angle<SoartingAngle) {
				int c_Ind = _wallFaces[wallFaceInd].layerCells[i];
				//tangential velocity in a cell
				double U_tan = U[c_Ind].rou/U[c_Ind].ro;
				if(U_tan > U_top[wallFaceInd])
				{
					U_top[wallFaceInd] = U_tan;
					CellSequence.push_back(c_Ind);
					v.value.push_back(U_tan);
					v.x.push_back(wInfo.distance);
				}else break;
			};
			v.size = v.x.size();
		};
		double h = v.x[v.size-1];

		//compute usefull functions
		function BM_FUN_t_visc = ComputeTurbViscosity(CellSequence, v);
		function BM_FUN_t_dif = ComputeTurbDiffusion(CellSequence, v);
		function lh_vel = ComputeLeeHarshaVelocity(v);

		BM_FUN_t_visc.WriteData("t_visc.dat");
		BM_FUN_t_dif.WriteData("t_def.dat");
		lh_vel.WriteData("t_LeeHarsha.dat");
		v.WriteData("t_velocity.dat");
		
		return;
	};
	function SolverTest(function& _t_vel)
	{
		//t_vel - function
		double U_inf = _t_vel.value[_t_vel.size - 1];
		double h = _t_vel.x[_t_vel.size - 1];
		BM_PAR_eps = medium.Viscosity/(h*U_inf);
		//BM_PAR_eps = 7.0e-4;
	
		BM_FUN_vel_der.value = BM_FUN_vel_der.value*(h/U_inf);
		BM_FUN_vel_der.x = BM_FUN_vel_der.x/h;

		BM_FUN_t_visc.value = BM_FUN_t_visc.value/(U_inf*h);
		BM_FUN_t_visc.x = BM_FUN_t_visc.x/h;

		BM_FUN_t_dif.value = BM_FUN_t_dif.value/(U_inf*h);
		BM_FUN_t_dif.x = BM_FUN_t_dif.x/h;

		BM_FUN_t_dif_der.value = BM_FUN_t_dif_der.value/(U_inf);
		BM_FUN_t_dif_der.x = BM_FUN_t_dif_der.x/h;

		BM_FUN_lh_vel_der.value = BM_FUN_lh_vel_der.value*(h/U_inf);
		BM_FUN_lh_vel_der.x = BM_FUN_lh_vel_der.x/h;

		BM_FUN_lh_vel_der.WriteData("BM_FUN_lh_vel_der.dat");
		BM_FUN_t_dif.WriteData("BM_FUN_t_dif.dat");
		BM_FUN_t_dif_der.WriteData("BM_FUN_t_dif_der.dat");
		BM_FUN_t_visc.WriteData("BM_FUN_t_visc.dat");
		_t_vel.WriteData("BM_FUN_t_vel.dat");
		BM_FUN_vel_der.WriteData("BM_FUN_vel_der.dat");

		return ComputeReinoldsStresses();
	};
	void TViscTest()
	{
		ComputeTurbStresses();
		SaveToTechPlot("t_viscocityBM.dat");
	};
};

#endif