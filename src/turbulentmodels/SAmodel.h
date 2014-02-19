#include "model.h"

struct SAConstansts {
	double k; //Karman constant
	double Cb1;
	double Cb2;
	double PrTurb; //Turbulent Prandtl number	
};


template<class RS> 
class SAModel : public Model<RS> {
protected:
	//Base model values
	DistributedEntityManager<RealValue> SAVariable;	//SA variable in each cell

	//Turbulent viscosity
	DistributedEntityManager<RealValue> TurbulentViscosityCells;	//Turbulent viscosity in each cell
	DistributedEntityManager<RealValue> TurbulentViscosityFaces;	//Turbulent viscosity on each face		

	//SA coefficients
	SAConstansts constants;	//Model constant parameters	

	//SA functions	

	//SA data structures
	std::map<int, double> _SAfluxes;		//Face.GlobalIndex -> convective flux
	std::map<int, double> _SAviscousFluxes;	//Face.GlobalIndex -> viscous flux
	std::map<int, double> _SAsourceTerms;	//Cell.GlobalIndex -> source term
	std::map<int, double> _boundarySAValues;//Modified turbulent viscosity value at boundaries	
	//Gradients
	std::map<int, Vector> _SAVariableGradientCells;	//SA variable gradient in cells
public:	
	//Constructor to assign default constant values
	SAModel() {
		constants.PrTurb = 0.677;
		constants.k = 0.4187;
		constants.Cb1 = 0.1355;
		constants.Cb2 = 0.622;

	};

	//Initialize data structures
	void Init() {		
		//Allocate memory for SAVariable
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for ( Cell* cell : cells) {
			SAVariable.add(RealValue(cell->GlobalIndex));
		};					
	};

	//Main steps of turbulence modeling
	void ComputeConvectiveSAFluxes() {
		std::vector<Face*> faces = _grid.faces.getLocalNodes();		

		//Compute convective flux for each cell face and apply boundary conditions								
		#pragma omp for
		for (Face* f : faces) {			
			RealValue SALeft = SAVariable[f->FaceCell_1];
			RealValue SARight;
			Vector velocityL = GetVelocity(U[f->FaceCell_1]);
			Vector velocityR;			
			if (f->isExternal) {
				SARight.value = GetBoundaryValues(*f);
				velocityR = GetVelocity(GetDummyCellValues(U[f->FaceCell_1], *f));
			} else {
				SARight = SAVariable[f->FaceCell_2];
				velocityR = GetVelocity(U[f->FaceCell_2]);
			};			

			//First order upwind			
			Vector velocityFace = 0.5 * (velocityL + velocityR);
			double velocity = f->FaceNormal * velocityFace;
			double flux = 0;
			if (velocity > 0) {
				flux = velocity * SALeft.value;
			} else {
				flux = velocity * SARight.value;
			};

			//Store fluxes
			_SAfluxes[f->GlobalIndex] = f->FaceSquare * flux;						
		};	

	};

	double GetSAVariable() {};
	void ComputeSAGradients() {
		//Compute gradients in cells
		ComputeFunctionGradient<RealValue>( &_SAVariableGradientCells, &SAVariable, GetSAVariable);
	};

	void ComputeViscousSAFluxes() {

	};


	//SA source terms
	double ProductionTerm() {
		return 0;
	};
	double DestructionTerm() {
		return 0;
	};
	void ComputeSASourceTerm() {

	};	

	//Boundary conditions
	void SetBoundarySAValue(std::string bcName, double value) {
		if (_grid.patchesNames.find(bcName) == _grid.patchesNames.end()) throw Exception("No such boundary availible");
		int bcMarker = _grid.patchesNames[bcName];
		_boundarySAValues[bcMarker] = value; 		
	};

	//Get boundary values depend on boundary conditions set
	double GetBoundaryValues(Face& face) {
		if (_boundarySAValues.find(face.BCMarker) != _boundarySAValues.end()) {
			return _boundarySAValues[face.BCMarker];
		} else {
			//If values not set assume zero turbulence level
			return 0;
		};
	};

	//Initial conditions
	void SetSAInitialConditions(double initSAValue) {
		std::vector<Cell*> cells = _grid.cells.getLocalNodes();
		for ( Cell& cell : cells) {
			SAVariable[cell->GlobalIndex] = initSAValue;
			SAVariable[cell->GlobalIndex].GlobalIndex = cell->GlobalIndex;
		};
		return;
	};

	//Main SA step
	void SAStep() {
		ComputeConvectiveSAFluxes();
		ComputeViscousSAFluxes();
		//ComputeSourceTerms();
		//DistributeFluxes();
	};


	//Load solution from CGNS file
	void ReadSolutionFromCGNS(std::string fname) {
		// TO DO process errors
		int fn;
		cg_open(fname.c_str(), CG_MODE_READ, &fn);

		// we know there is only one base (real working code would check!)
		int index_base=1;

		/* Check the cell and physical dimensions of the base. */	
		int physDim;
		int cellDim;
		char cgnsName[255];
		if(cg_base_read(fn, index_base, cgnsName, &cellDim, &physDim) != CG_OK) cg_error_exit();
		
		// we know there is only one zone (real working code would check!)
		int index_zone=1;
		// we know there is only one FlowSolution_t (real working code would check!)
		int index_flow=1;
		// get zone size (and name - although not needed here)
		cgsize_t sizes[3];
		char zone_name[255];
		cg_zone_read(fn, index_base, index_zone, zone_name, sizes);

		// check if size is consistent with grid size
		if (sizes[1] != U.size()) {
			std::cout <<"Error loading solution, wrong number of solution records!\n";
			return;
		};

		//   upper range index - use vertex dimensions
		//  checking GridLocation first (real working code would check
		//  to make sure there are no Rind cells also!):
		CGNS_ENUMT(GridLocation_t) location;
		char solution_name[255];
		cg_sol_info(fn,index_base,index_zone,index_flow, solution_name, &location);
		if (location != CellCenter) {
			std::cout <<"Error loading solution, GridLocation must be CellCenter!\n";
			return;
		};  
		// lower range index
		cgsize_t irmin[1];
		irmin[0] = 1;
		cgsize_t irmax[1];
		irmax[0] = U.size();
		// read flow solution
		// allocate input buffer for Double precision values
		double* inputBuffer = new double[U.size()];

		// load density
		cg_field_read(fn,index_base,index_zone,index_flow, "Density" , RealDouble, irmin, irmax, (void*)inputBuffer);
		for (int i = 0; i<U.size(); i++) {
			U[i+1].ro = inputBuffer[i];
		};

		// load velocity components
		cg_field_read(fn,index_base,index_zone,index_flow, "VelocityX" , RealDouble, irmin, irmax, (void*)inputBuffer);
		for (int i = 0; i<U.size(); i++) {
			U[i+1].rou = U[i+1].ro * inputBuffer[i];
		};

		// load Y component in case its not 1D problem
		for (int i = 0; i<U.size(); i++) U[i+1].rov = 0;			
		if (physDim > 1) {
			cg_field_read(fn,index_base,index_zone,index_flow, "VelocityY" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].rov = U[i+1].ro * inputBuffer[i];
			};
		};

		// load Z component in case its not 1D or 2D problem
		for (int i = 0; i<U.size(); i++) U[i+1].row = 0;			
		if (physDim > 2) {
			cg_field_read(fn,index_base,index_zone,index_flow, "VelocityZ" , RealDouble, irmin, irmax, (void*)inputBuffer);
			for (int i = 0; i<U.size(); i++) {
				U[i+1].row = U[i+1].ro * inputBuffer[i];
			};
		};

		// load stagnation energy		
		cg_field_read(fn,index_base,index_zone,index_flow, "EnergyStagnation" , RealDouble, irmin, irmax, (void*)inputBuffer);  
		for (int i = 0; i<U.size(); i++) {
			U[i+1].roE = U[i+1].ro * inputBuffer[i];
		};
		
		// load SA variable
		cg_field_read(fn,index_base,index_zone,index_flow, "Modified_Turbulent_Viscosity" , RealDouble, irmin, irmax, (void*)inputBuffer);  
		for (int i = 0; i<U.size(); i++) {
			SAVariable[i+1].value = inputBuffer[i];
		};
		

		// free buffer memory
		delete[] inputBuffer;

		// close CGNS file		
		cg_close(fn);
	};

	//Write solution to CGNS file
	void WriteSolutionToCGNS(std::string fname) {

	};

	//Postprocessing
	

//End of class SAModel
};