#include "model.h"

struct CellVariable
{
	double u;
	double v;

	CellVariable(double _u, double _v):u(_u), v(_v) {};
	CellVariable(){};
};
struct Incompress2DFlowProperties{
	double U_inf;
	double density;
	double Temperature;
	double Gamma;
	double Cv;
	double viscosity;	//dinamic viscosity

	void Set_Uinf(double _U_inf) {
		U_inf = _U_inf;
	};
	void Set_density(double _density)
	{
		density = _density;
	};
	void Set_Temperature(double _Temperature)
	{
		Temperature = _Temperature;
	}
	void Set_Gamma(double _Gamma)
	{
		Gamma = _Gamma;
	};
	void Set_Cv(double _Cv)
	{
		Cv = _Cv;
	};
	void Set_viscosity(double _viscosity)
	{
		viscosity = _viscosity;
	};
};
////functions and parameters for Levenberg–Marquardt solver////
//compute residual for V velocity and write it in f array
struct Flow2DTaskParams
{
	Grid* g;
	int N,M;
	double L_x,L_y;
	Incompress2DFlowProperties medium;
	std::vector<std::vector<int>> CList;
	std::map<int, CellVariable> U;
};

//residual of levenberg-markquardt
function Flow2DRes;

void Flow2DResWrite(char* filename) {
	for(int i=0; i<Flow2DRes.value.size(); i++)	Flow2DRes.x.push_back(i);
	Flow2DRes.WriteData(filename);
	return;
};


void ComputeWriteResV(const alglib::real_1d_array &x, alglib::real_1d_array &f, void *ptr)
{
	Flow2DTaskParams* taskParams = (Flow2DTaskParams*)ptr;
	Grid* g = taskParams->g;
	int N = taskParams->N;
	int M = taskParams->M;
	Incompress2DFlowProperties medium = taskParams->medium;
	double res = 0;

	for(int i=0; i<N; i++)
	{
		//cells sizes
		double h_x_l, h_x_r, h_y_u, h_y_d;;
		CellVariable Ul, Ur, Uu, Ud;

		for(int j=0; j<M; j++)
		{
			CellVariable U;
			U.u = taskParams->U[taskParams->CList[i][j]].u;
			U.v = x[taskParams->CList[i][j]];
			//bottom
			if(j==0) {
				h_y_d = 2.0 * g->cells[taskParams->CList[i][j]].CellCenter.y;
				Ud = CellVariable(-U.u, -U.v);
			}else {
				h_y_d = g->cells[taskParams->CList[i][j]].CellCenter.y - g->cells[taskParams->CList[i][j-1]].CellCenter.y;
				Ud.u = taskParams->U[taskParams->CList[i][j-1]].u;
				Ud.v = x[taskParams->CList[i][j-1]];
			};
			//top
			if(j==M-1) {
				h_y_u = 2.0 * (taskParams->L_y - g->cells[taskParams->CList[i][j]].CellCenter.y);
				Uu = CellVariable(medium.U_inf, U.v);
			}else {
				h_y_u = g->cells[taskParams->CList[i][j+1]].CellCenter.y - g->cells[taskParams->CList[i][j]].CellCenter.y;
				Uu.u = taskParams->U[taskParams->CList[i][j+1]].u;
				Uu.v = x[taskParams->CList[i][j+1]];
			};
			double hy = h_y_u + h_y_d;
			//left
			if(i==0) {
				h_x_l = 2.0 * g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ul = CellVariable(medium.U_inf, 0);
			}else {
				h_x_l = g->cells[taskParams->CList[i][0]].CellCenter.x - g->cells[taskParams->CList[i-1][0]].CellCenter.x;
				Ul.u = taskParams->U[taskParams->CList[i-1][j]].u;
				Ul.v = x[taskParams->CList[i-1][j]];
			};
			//right
			if(i==N-1) {
				h_x_r = 2.0 * (taskParams->L_x - g->cells[taskParams->CList[i][0]].CellCenter.x);
				Ur = U;
			}else {
				h_x_r = g->cells[taskParams->CList[i+1][0]].CellCenter.x - g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ur.u = taskParams->U[taskParams->CList[i+1][j]].u;
				Ur.v = x[taskParams->CList[i+1][j]];
			};
			double hx = h_x_l + h_x_r;

			double dUdx = (Ur.u - Ul.u)/hx;
			double dVdy = (Uu.v - Ud.v)/hy;

			//compute a components of optimized function
			int k = i*M + j;
			f[k] = dUdx + dVdy;
			res += f[k]*f[k];
		};
	};
	Flow2DRes.value.push_back(res);
	return; 
};
void ComputeWriteResU(const alglib::real_1d_array &x, alglib::real_1d_array &f, void *ptr)
{
	Flow2DTaskParams* taskParams = (Flow2DTaskParams*)ptr;
	Grid* g = taskParams->g;
	int N = taskParams->N;
	int M = taskParams->M;
	Incompress2DFlowProperties medium = taskParams->medium;
	double res = 0;

	for(int i=0; i<N; i++)
	{
		//cells sizes
		double h_x_l, h_x_r, h_y_u, h_y_d;;
		CellVariable Ul, Ur, Uu, Ud;

		for(int j=0; j<M; j++)
		{
			CellVariable U;
			U.u = x[taskParams->CList[i][j]];
			U.v = taskParams->U[taskParams->CList[i][j]].v;
			//bottom
			if(j==0) {
				h_y_d = 2.0 * g->cells[taskParams->CList[i][j]].CellCenter.y;
				Ud = CellVariable(-U.u, -U.v);
			}else {
				h_y_d = g->cells[taskParams->CList[i][j]].CellCenter.y - g->cells[taskParams->CList[i][j-1]].CellCenter.y;
				Ud.u = x[taskParams->CList[i][j-1]];
				Ud.v = taskParams->U[taskParams->CList[i][j-1]].v;
			};
			//top
			if(j==M-1) {
				h_y_u = 2.0 * (taskParams->L_y - g->cells[taskParams->CList[i][j]].CellCenter.y);
				Uu = CellVariable(medium.U_inf, U.v);
			}else {
				h_y_u = g->cells[taskParams->CList[i][j+1]].CellCenter.y - g->cells[taskParams->CList[i][j]].CellCenter.y;
				Uu.u = x[taskParams->CList[i][j+1]];
				Uu.v = taskParams->U[taskParams->CList[i][j+1]].v;
			};
			double hy = h_y_u + h_y_d;
			//left
			if(i==0) {
				h_x_l = 2.0 * g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ul = CellVariable(medium.U_inf, 0);
			}else {
				h_x_l = g->cells[taskParams->CList[i][0]].CellCenter.x - g->cells[taskParams->CList[i-1][0]].CellCenter.x;
				Ul.u = x[taskParams->CList[i-1][j]];
				Ul.v = taskParams->U[taskParams->CList[i-1][j]].v;
			};
			//right
			if(i==N-1) {
				h_x_r = 2.0 * (taskParams->L_x - g->cells[taskParams->CList[i][0]].CellCenter.x);
				Ur = U;
			}else {
				h_x_r = g->cells[taskParams->CList[i+1][0]].CellCenter.x - g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ur.u = x[taskParams->CList[i+1][j]];
				Ur.v = taskParams->U[taskParams->CList[i+1][j]].v;
			};
			double hx = h_x_l + h_x_r;

			double dUdx = (Ur.u - Ul.u)/hx;
			double dUdy = (Uu.u - Ud.u)/hy;
			double dUdy2 = 4*(Uu.u + Ud.u - 2*U.u)/(hy*hy);

			//compute a components of optimized function
			int k = i*M + j;
			double ny = medium.viscosity/medium.density;
			f[k] = U.u*dUdx + U.v*dUdy - ny*dUdy2;
			res += f[k]*f[k];
		};
	};
	Flow2DRes.value.push_back(res);
	return; 
};
void ComputeWriteResB(const alglib::real_1d_array &x, alglib::real_1d_array &f, void *ptr)
{
	Flow2DTaskParams* taskParams = (Flow2DTaskParams*)ptr;
	Grid* g = taskParams->g;
	int N = taskParams->N;
	int M = taskParams->M;
	int Nc = N*M;		//cells number
	Incompress2DFlowProperties medium = taskParams->medium;
	double res = 0;

	for(int i=0; i<N; i++)
	{
		//cells sizes
		double h_x_l, h_x_r, h_y_u, h_y_d;;
		CellVariable Ul, Ur, Uu, Ud;

		for(int j=0; j<M; j++)
		{
			CellVariable U;
			U.u = x[taskParams->CList[i][j]];
			U.v = x[taskParams->CList[i][j] + Nc];
			//bottom
			if(j==0) {
				h_y_d = 2*g->cells[taskParams->CList[i][j]].CellCenter.y;
				Ud = CellVariable(-U.u, -U.v);
				//Ud = CellVariable(0, 0);
			}else {
				h_y_d = g->cells[taskParams->CList[i][j]].CellCenter.y - g->cells[taskParams->CList[i][j-1]].CellCenter.y;
				Ud.u = x[taskParams->CList[i][j-1]];
				Ud.v = x[taskParams->CList[i][j-1] + Nc];
			};
			//top
			if(j==M-1) {
				h_y_u = 2.0 * (taskParams->L_y - g->cells[taskParams->CList[i][j]].CellCenter.y);
				Uu = CellVariable(medium.U_inf, U.v);
			}else {
				h_y_u = g->cells[taskParams->CList[i][j+1]].CellCenter.y - g->cells[taskParams->CList[i][j]].CellCenter.y;
				Uu.u = x[taskParams->CList[i][j+1]];
				Uu.v = x[taskParams->CList[i][j+1] + Nc];
			};
			double hy = h_y_u + h_y_d;
			//left
			if(i==0) {
				h_x_l = 2.0 * g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ul = CellVariable(medium.U_inf, 0);
			}else {
				h_x_l = g->cells[taskParams->CList[i][0]].CellCenter.x - g->cells[taskParams->CList[i-1][0]].CellCenter.x;
				Ul.u = x[taskParams->CList[i-1][j]];
				Ul.v = x[taskParams->CList[i-1][j] + Nc];
			};
			//right
			if(i==N-1) {
				h_x_r = 2.0 * (taskParams->L_x - g->cells[taskParams->CList[i][0]].CellCenter.x);
				Ur = U;
			}else {
				h_x_r = g->cells[taskParams->CList[i+1][0]].CellCenter.x - g->cells[taskParams->CList[i][0]].CellCenter.x;
				Ur.u = x[taskParams->CList[i+1][j]];
				Ur.v = x[taskParams->CList[i+1][j] + Nc];
			};
			double hx = h_x_l + h_x_r;

			//for uniform grid
			//double dUdx = (Ur.u - Ul.u)/hx;
			//double dUdy = (Uu.u - Ud.u)/hy;
			//double dUdy2 = 4*(Uu.u + Ud.u - 2*U.u)/(hy*hy);
			//double dVdy = (Uu.v - Ud.v)/hy;
			//for not uniform grid
			double dUdx = (h_x_r - h_x_l)*U.u/(h_x_r*h_x_l);
			dUdx += h_x_l*Ur.u/(h_x_r*hx) - h_x_r*Ul.u/(h_x_l*hx);
			double dUdy = (h_y_u - h_y_d)*U.u/(h_y_u*h_y_d);
			dUdy += h_y_d*Uu.u/(h_y_u*hy) - h_y_u*Ud.u/(h_y_d*hy);
			double dVdy = (h_y_u - h_y_d)*U.v/(h_y_u*h_y_d);
			dUdy += h_y_d*Uu.v/(h_y_u*hy) - h_y_u*Ud.v/(h_y_d*hy);
			double dUdy2 = Uu.u/(h_y_u*hy) + Ud.u/(h_y_d*hy) - U.u/(h_y_u*h_y_d);
			dUdy2 *= 2.0;

			//compute a components of optimized function
			int k1 = 2*(i*M + j);
			int k2 = k1+1;
			double ny = medium.viscosity/medium.density;
			f[k1] = U.u*dUdx + U.v*dUdy - ny*dUdy2;
			f[k2] = dUdx + dVdy;
			res += f[k1]*f[k1] + f[k2]*f[k2];
		};
	};
	Flow2DRes.value.push_back(res);
	return; 
};

class Flow2D
{
public:
	Grid _grid;	//2D grid
	Incompress2DFlowProperties medium;	//flow and medium parameters
	double kinetic_visc;	//kinetic_viscosity
	std::map<int, CellVariable> U;	//from cell to variables map
	int N_x, M_y;	//x and y direction cells numbers
	double L_x, L_y;	//x and y length

	//constructors
	Flow2D(){}; 
	Flow2D(int _N_x, int _M_y):N_x(_N_x), M_y(_M_y){};

	//make convectivity list
	std::vector<std::vector<int>> CList;
	void CreateConvectivityList()
	{
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		CList.resize(N_x);	
		//main cycle
		for(int i=0; i<cells.size();) {
			for(int j=0; j<N_x; j++)	{
				CList[j].push_back(cells[i]->GlobalIndex);
				i++;
			};
		};
	};
	
	//bind U to grid
	void BindGrid(Grid& grid)	
	{
		U.clear();
		_grid = grid;
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++)
		{
			Cell c = *cells[i];
			U[c.GlobalIndex] = CellVariable(0,0);
		};
		return;
	};
	//Set M and N
	void SetGridParameters(int N, int M, double Lx, double Ly)
	{
		N_x = N;
		M_y = M;
		L_x = Lx;
		L_y = Ly;
	};
	//compute kinetic_visc
	void ComputeKineticViscosity()
	{
		kinetic_visc = medium.viscosity/medium.density;
	};
	//set initial condition close to Blasius
	void SetInitialCondition(double x_plate)
	{
		CellVariable FreeStream(medium.U_inf, 0);

		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for(int i=0; i<cells.size(); i++) {
			Vector R = cells[i]->CellCenter;
			if(R.x<x_plate) {
				U[cells[i]->GlobalIndex] = FreeStream;
				continue;
			}else {
			double nu = R.y*sqrt(medium.U_inf/((R.x - x_plate)*kinetic_visc));
			if(nu>3.2)  {
				U[cells[i]->GlobalIndex] = FreeStream;
				continue;
			};
			U[cells[i]->GlobalIndex].u = nu*medium.U_inf/3.2;
			U[cells[i]->GlobalIndex].v = 0;
			};
		};
	};

	//solver functions
	//compute V velocity field from conservation of mass
	void ComputeVfromU() {
		Flow2DTaskParams TP;	//to pass parameters to target functions
		TP.g = &_grid;
		TP.CList = CList;
		TP.M = M_y; TP.N = N_x;
		TP.medium = medium;
		TP.U = U;
		TP.L_x = L_x;
		TP.L_y = L_y;

		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		alglib::real_1d_array v;	//initial approach
		v.setlength(cells.size());
		for(int i=0; i<v.length(); i++)
		{
			v[cells[i]->GlobalIndex] = U[cells[i]->GlobalIndex].v;
		};
				
		double epsg = 0.00001;
		double epsf = 0;
		double epsx = 0;
		alglib::ae_int_t maxits = 0;
		alglib::ae_int_t size = v.length();
		alglib::minlmstate state;
		alglib::minlmreport rep;
		double diffstep = 0.0001;

		alglib::minlmcreatev(size, v, diffstep, state);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
		alglib::minlmoptimize(state, ComputeWriteResV, NULL, (void *)&TP);
		alglib::minlmresults(state, v, rep);

		//write result
		for(int i=0; i<v.length(); i++) U[cells[i]->GlobalIndex].v = v[i];
		return;
	};
	void ComputeUfromV() {
		Flow2DTaskParams TP;	//to pass parameters to target functions
		TP.g = &_grid;
		TP.CList = CList;
		TP.M = M_y; TP.N = N_x;
		TP.medium = medium;
		TP.U = U;
		TP.L_x = L_x;
		TP.L_y = L_y;

		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		alglib::real_1d_array u;	
		u.setlength(cells.size());
		for(int i=0; i<u.length(); i++)
		{
			u[cells[i]->GlobalIndex] = U[cells[i]->GlobalIndex].u;	//initial approach
		};
				
		double epsg = 0.0001;
		double epsf = 0;
		double epsx = 0;
		alglib::ae_int_t maxits = 0;
		alglib::ae_int_t size = u.length();
		alglib::minlmstate state;
		alglib::minlmreport rep;
		double diffstep = 0.0001;

		alglib::minlmcreatev(size, u, diffstep, state);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
		alglib::minlmoptimize(state, ComputeWriteResU, NULL, (void *)&TP);
		alglib::minlmresults(state, u, rep);
		std::cout << rep.terminationtype << '\n';

		//write result
		for(int i=0; i<u.length(); i++) U[cells[i]->GlobalIndex].u = u[i];
		return;
	};
	void ComputeBlasius() {
		Flow2DTaskParams TP;	//to pass parameters to target functions
		TP.g = &_grid;
		TP.CList = CList;
		TP.M = M_y; TP.N = N_x;
		TP.medium = medium;
		TP.U = U;
		TP.L_x = L_x;
		TP.L_y = L_y;

		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		alglib::real_1d_array x;	
		x.setlength(2*cells.size());
		for(int i=0; i<cells.size(); i++)
		{
			x[cells[i]->GlobalIndex] = U[cells[i]->GlobalIndex].u;	//initial approach
		};
		for(int i=0; i<cells.size(); i++)
		{
			x[cells[i]->GlobalIndex + cells.size()] = U[cells[i]->GlobalIndex].v;	//initial approach
		};			
				
		double epsg = 0.1;
		double epsf = 0;
		double epsx = 0;
		alglib::ae_int_t maxits = 200;
		alglib::ae_int_t size = x.length();
		alglib::minlmstate state;
		alglib::minlmreport rep;
		double diffstep = 0.0001;

		alglib::minlmcreatev(size, x, diffstep, state);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
		alglib::minlmoptimize(state, ComputeWriteResB, NULL, (void *)&TP);
		alglib::minlmresults(state, x, rep);
		std::cout << rep.terminationtype << '\n';

		//write result
		for(int i=0; i<cells.size(); i++) {
			U[cells[i]->GlobalIndex].u = x[i];
		};
		for(int i=0; i<cells.size(); i++) {
			U[cells[i]->GlobalIndex].v = x[i + cells.size()];
		};
		return;
	};

	//Save solution to TecPlot
	void SaveToTechPlot(std::string fname) {
		std::ofstream ofs(fname);
		ofs<<std::scientific;

		//Header
		ofs<<"VARIABLES= \"X\", \"Y\", \"Rho\", \"u\", \"v\", \"T\", \"P\"";
		ofs<<"\n";
			
		ofs<<"ZONE T=\"D\"\n";
		ofs<<"N=" << _grid.nodes.size() << ", E=" << _grid.cells.size() <<", F=FEBLOCK, ET=QUADRILATERAL\n";

		ofs<<"VARLOCATION = (NODAL, NODAL, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED, CELLCENTERED";
		ofs<<")\n";

		//Map all node global indexes to natural numbers
		std::map<int,int> toNaturalIndex;
		std::set<int> nodeIndexes = _grid.nodes.getAllIndexes();
		int counter = 1;
		for (std::set<int>::iterator it = nodeIndexes.begin(); it != nodeIndexes.end(); it++) toNaturalIndex[*it] = counter++;
			
		//Access local nodes, faces, cells and flow data
		std::vector<Node*> nodes = _grid.nodes.getLocalNodes();
		std::vector<Face*> faces = _grid.faces.getLocalNodes();
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();

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
			ofs<< medium.density <<"\n";
		};
				
		//u
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<< U[idx].u <<"\n";
		};

		//v
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<< U[idx].v <<"\n";		
		};

		//T
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			ofs<<medium.Temperature<<"\n";
		};
		//P
		for (int i = 0; i<cells.size(); i++) {
			int idx = cells[i]->GlobalIndex;
			double P = (medium.Gamma-1.0) * medium.density * medium.Gamma * medium.Cv * medium.Temperature;
			ofs<<P<<"\n";
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
	//Save solution in a file //TO DO not portable
    void SaveSolution(std::string fname = "SavedSolution.sol")
        {
                std::ofstream ofs;                
                ofs.open(fname, std::ios::binary | std::ios::out);                
                std::vector<Cell*> cells = _grid.cells.getLocalNodes();                                
				ofs.write( reinterpret_cast<char*>( &medium.density ), sizeof medium.density);
                for(int i=0; i<cells.size(); i++)
                {
                        int Index = cells[i]->GlobalIndex;        //Global Index of cell
                        ofs.write( reinterpret_cast<char*>( &Index ), sizeof Index );
						ofs.write( reinterpret_cast<char*>( &U[Index].u ), sizeof U[Index].u );
                        ofs.write( reinterpret_cast<char*>( &U[Index].v ), sizeof U[Index].v );                        
                };
                ofs.close();
        };
    //Load a solution file //TO DO not portable
    void LoadSolution(std::string fname = "SavedSolution.txt")
        {
                std::ifstream ifs;
                ifs.open(fname, std::ios::binary | std::ios::in);
				ifs.read( reinterpret_cast<char*>( &medium.density ), sizeof medium.density);
                while(!ifs.eof())
                {
                        int Index;                        
                        ifs.read( reinterpret_cast<char*>( &Index ), sizeof Index );
                        ifs.read( reinterpret_cast<char*>( &U[Index].u ), sizeof U[Index].u );
                        ifs.read( reinterpret_cast<char*>( &U[Index].v ), sizeof U[Index].v );                       
                };
                ifs.close();
        };
};


/*
	//Patches and corresponding functions
	class BoundaryCondition {		
	public:
		Flow2D& Fl2D;
				
		BoundaryCondition(Flow2D& _Fl2D) : Fl2D(_Fl2D) {};
		virtual CellVariable getDummyValues(CellVariable UL, const Face& face) = 0;		
	};
	class NoSlipBoundaryCondition : public BoundaryCondition {
	public:
		NoSlipBoundaryCondition(Flow2D& _Fl2D):BoundaryCondition(_Fl2D){};
		CellVariable getDummyValues(CellVariable UL, const Face& face) {
			UL.u *= -1;
			UL.v *= -1;
			return UL;
		};
	};
	class InletBoundaryCondition : public BoundaryCondition {
	public:
		InletBoundaryCondition(Flow2D& _Fl2D):BoundaryCondition(_Fl2D){};
		CellVariable getDummyValues(CellVariable UL, const Face& face) {
			UL.u = Fl2D.medium.U_inf;
			UL.v *= 0;
			return UL;
		};
	};
	class OutletBoundaryCondition : public BoundaryCondition {
	public:
		OutletBoundaryCondition(Flow2D& _Fl2D):BoundaryCondition(_Fl2D){};
		CellVariable getDummyValues(CellVariable UL, const Face& face) {
			return UL;
		};
	};
	class ConstantVelocityXBoundaryCondition : public BoundaryCondition {
	public:
		ConstantVelocityXBoundaryCondition(Flow2D& _Fl2D):BoundaryCondition(_Fl2D){};
		CellVariable getDummyValues(CellVariable UL, const Face& face) {
			UL.u = Fl2D.medium.U_inf;
			return UL;
		};
	};
	class SymmetryXBoundaryCondition : public BoundaryCondition {
	public:
		SymmetryXBoundaryCondition(Flow2D& _Fl2D):BoundaryCondition(_Fl2D){};
		CellVariable getDummyValues(CellVariable UL, const Face& face) {
			UL.u = Fl2D.medium.U_inf;
			UL.v *= -1.0;
			return UL;
		};
	};
	void SetBoundaryCondition(std::string BCName, BoundaryCondition& bc) {
		if (_grid.patchesNames.find(BCName) == _grid.patchesNames.end()) throw Exception("No such boundary availible");
		int bcMarker = _grid.patchesNames[BCName];
		_boundaryConditions[bcMarker] = &bc; 
	};
	//Boundary conditions map
	std::map<int, BoundaryCondition*> _boundaryConditions;
*/