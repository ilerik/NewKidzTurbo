#ifndef TURBO_MAINSOLVER_PARALLELHELPER
#define TURBO_MAINSOLVER_PARALLELHELPER

#include <mpi.h>
#include <time.h>
#include <assert.h>

//Class that implements MPI-based operations for distributed memory architecture
class ParallelHelper {
	bool isInitilized; //was MPI initialization successful
	int _nProcessors; //total number of processes
	int _rank;  //current processor rank in _comm
	MPI_Comm _comm;	//global communicator
public:
	//Accessor functions
	inline MPI_Comm getComm() {
		return _comm;
	};

	inline int getRank() {
		return _rank;
	};

	inline int getProcessorNumber() {
		return _nProcessors;
	};

	//Constructor
	ParallelHelper() {
		isInitilized = false;
	};

	//Initialize MPI programm
	void Init(int *argc, char **argv[]) {
		_comm = MPI_COMM_WORLD;
		MPI_Init(argc, argv);		
		MPI_Comm_size(_comm, &_nProcessors);
		MPI_Comm_rank(_comm, &_rank);
		isInitilized = true;
	};

	//Finilize MPI programm
	void Finalize() {
		MPI_Finalize();
		isInitilized = false;
	};

	//Is master node
	inline bool IsMaster() {
		return _rank == 0;
	};

	//Barrier
	double Barrier() {
		double t = 0.0;
		clock_t start, stop;
		/* Start timer */
		assert((start = clock())!=-1);
		MPI_Barrier(_comm);

		/* Stop timer */
		stop = clock();
		t = (double) (stop-start)/CLOCKS_PER_SEC;
		return t;		
	};

	//Gather one integer from each process
	void GatherCounts(int N, std::vector<int>& result) {		
		if (IsMaster()) {			
			result.resize(_nProcessors);
			MPI_Gather(&N, 1, MPI_INT, &result[0], 1, MPI_INT, 0, _comm);
		} else {
			MPI_Gather(&N, 1, MPI_INT, NULL, 1, MPI_INT, 0, _comm);
		};		
	};
	
	//All gather one integer
	void AllgatherCounts(int N, std::vector<int>& result) {				
		result.resize(_nProcessors);
		MPI_Allgather(&N, 1, MPI_INT, &result[0], 1, MPI_INT, _comm);			
	};

	//Gather arrays of double of different size on master node
	void GathervInt(std::vector<int>& local, std::vector<int>& counts, std::vector<int>& result) {		
		if (IsMaster()) {			
			//Make displacement array
			std::vector<int> displs(counts.size());
			int totalSize = 0;
			for (int i = 0; i<counts.size(); i++) {
				displs[i] = totalSize;
				totalSize += counts[i];
			};
			result.resize(totalSize);
			if (local.size() != 0) {
				MPI_Gatherv(&local[0], local.size(), MPI_INT, &result[0], &counts[0], &displs[0], MPI_INT, 0, _comm);
			} else {
				MPI_Gatherv(NULL, local.size(), MPI_INT, &result[0], &counts[0], &displs[0], MPI_INT, 0, _comm);
			};
		} else {
			if (local.size() != 0) {
				MPI_Gatherv(&local[0], local.size(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, _comm);
			} else {
				MPI_Gatherv(NULL, local.size(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, _comm);
			};			
		};		
	};

	//Gather arrays of double of different size on master node
	void GathervDouble(std::vector<double>& local, std::vector<int>& counts, std::vector<double>& result) {		
		if (IsMaster()) {			
			//Make displacement array
			std::vector<int> displs(counts.size());
			int totalSize = 0;
			for (int i = 0; i<counts.size(); i++) {
				displs[i] = totalSize;
				totalSize += counts[i];
			};
			result.resize(totalSize);
			if (local.size() != 0) {
				MPI_Gatherv(&local[0], local.size(), MPI_LONG_DOUBLE, &result[0], &counts[0], &displs[0], MPI_LONG_DOUBLE, 0, _comm);
			} else {
				MPI_Gatherv(NULL, local.size(), MPI_LONG_DOUBLE, &result[0], &counts[0], &displs[0], MPI_LONG_DOUBLE, 0, _comm);
			};			
		} else {
			if (local.size() != 0) {
				MPI_Gatherv(&local[0], local.size(), MPI_LONG_DOUBLE, NULL, NULL, NULL, MPI_LONG_DOUBLE, 0, _comm);
			} else {
				MPI_Gatherv(NULL, local.size(), MPI_LONG_DOUBLE, NULL, NULL, NULL, MPI_LONG_DOUBLE, 0, _comm);
			};					
		};		
	};

	//Wrapper for Allgatherv (int)
	void Allgatherv(std::vector<int>& local, std::vector<int>& counts, std::vector<int>& result) {
		//Make displacement array
		std::vector<int> displs(counts.size());
		int totalSize = 0;
		for (int i = 0; i<counts.size(); i++) {
			displs[i] = totalSize;
			totalSize += counts[i];
		};
		result.resize(totalSize);
		MPI_Allgatherv(&local[0], local.size(), MPI_INT, &result[0], &counts[0], &displs[0], MPI_INT, _comm);		
	};
	////Agregation
	
	//Sum of integers
	inline int SumInt(int x) {
		int res;
		MPI_Allreduce(&x, &res, 1, MPI_INT, MPI_SUM, _comm);
		return res;
	};

	//Sum of integers
	inline double SumDouble(double x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_SUM, _comm);
		return res;
	};

	//Minimum
	inline double Min(double x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_MIN, _comm);
		return res;
	};


public:
	//Problem oriented routines and data
	int nVariables; //Number of conservative variables
	std::vector<int> part; //Partitioning
	std::map<int, std::vector<double>> RequestedValues; //Requested values for cells
	std::map<int, std::vector<Vector>> RequestedGradients; //Requested gradients for cells

	//To send by processor
	std::vector<std::vector<int>> toSendValuesByProc;
	std::vector<int> toSendValuesNumberByProc;

	//To recieve by processor
	std::set<int> toRecvValues;
	std::vector<std::vector<int>> toRecvValuesByProc;
	std::vector<int> toRecvValuesNumberByProc;		

	//Set cells partitioning
	void SetCellsPartitioning(const std::vector<int>& part_) {
		part = part_;
	};

	//Clear all requests
	void ClearRequests() {
		//Clear values requests
		toSendValuesByProc.clear();
		toRecvValuesByProc.clear();
		toRecvValues.clear();

		//Resize values structures		
		toSendValuesByProc.resize(_nProcessors, std::vector<int>());
		toRecvValuesByProc.resize(_nProcessors, std::vector<int>());
	};

	//Add request
	void RequestValues(int globalIndex) {
		int processor = part[globalIndex];
		if (_rank == processor) return; //Skip local cell requests
		toRecvValues.insert(globalIndex);
	};

	//Initialize data structures for exchange
	void InitExchange() {				
		//Cell indexes per processor to recieve		
		for (int cellID : toRecvValues) {
			toRecvValuesByProc[part[cellID]].push_back(cellID);
		};

		//Determine number of needed values by each processor		
		toRecvValuesNumberByProc.resize(_nProcessors);
		for (int rank = 0; rank < _nProcessors; rank++) toRecvValuesNumberByProc[rank] = toRecvValuesByProc[rank].size();
		toSendValuesNumberByProc.resize(_nProcessors);

		//Exchange amounts of needed cells		
		MPI_Alltoall(&toRecvValuesNumberByProc[0], 1, MPI_INT, &toSendValuesNumberByProc[0], 1, MPI_INT, _comm);

		//Echange indexes of cells to send
		int s = 0;
		int r = 0;
		std::vector<int> sdispl;
		std::vector<int> rdispl;
		std::vector<int> recvbuf;
		std::vector<int> sendbuf;		
		for (int i = 0; i<_nProcessors; i++) {
			sdispl.push_back(s);
			rdispl.push_back(r);
			//From current proc to i-th proc send what i want to recieve			
			for (int j = 0; j<toRecvValuesNumberByProc[i]; j++) {
				sendbuf.push_back(toRecvValuesByProc[i][j]);
			};
			s += toRecvValuesNumberByProc[i];
			r += toSendValuesNumberByProc[i];			
		};
		//Make sure that sendbuf isnt empty
		sendbuf.push_back(0);

		//Allocate memory to recieve what cells j-th proc want me to send
		recvbuf.resize(r + 1);

		//Exchange cell indexes		
		MPI_Alltoallv(&sendbuf[0], &toRecvValuesNumberByProc[0], &sdispl[0], MPI_INT, &recvbuf[0], &toSendValuesNumberByProc[0], &rdispl[0], MPI_INT, _comm);

		//Extract cell indexes to send
		for (int i = 0; i < _nProcessors; i++) {
			for (int j = rdispl[i]; j < rdispl[i] + toSendValuesNumberByProc[i]; j++) {
				toSendValuesByProc[i].push_back(recvbuf[j]);
			};
		};
	};

	//Exchange values
	void ExchangeValues(Grid& grid, std::vector<double>& values) {
		printf("00");
		//Allocate memory and fill datastructures
		int s = 0;
		int r = 0;
		std::vector<int> sdispl;
		std::vector<int> rdispl;
		std::vector<int> recvbuf;
		std::vector<int> sendbuf;		
		for (int i = 0; i<_nProcessors; i++) {
			sdispl.push_back(s);
			rdispl.push_back(r);
			//From current proc to i-th proc send values that proc requested			
			for (int j = 0; j<toSendValuesNumberByProc[i]; j++) {
				//Cell index
				int cellGlobalIndex = toSendValuesByProc[i][j];
				int cellLocalIndex = grid.cellsGlobalToLocal[cellGlobalIndex];
				//Write values from cell to send buffer
				for (int k = 0; k < nVariables; k++) {
					sendbuf.push_back(values[cellLocalIndex * nVariables + k]);
				};
			};
			s += toSendValuesNumberByProc[i] * nVariables;
			r += toRecvValuesNumberByProc[i] * nVariables;			
		};

		//Allocate recieve buffer
		recvbuf.resize(r);
		printf("11");

		//Exchange values
		MPI_Alltoallv(&sendbuf[0], &toSendValuesNumberByProc[0], &sdispl[0], MPI_LONG_DOUBLE, &recvbuf[0], &toRecvValuesNumberByProc[0], &rdispl[0], MPI_LONG_DOUBLE, _comm);		

		//Write recieved values to appropriate data structure
		int pointer = 0;
		for (int i = 0; i<_nProcessors; i++) {								
			for (int j = 0; j<toRecvValuesNumberByProc[i]; j++) {			
				//Cell index
				int cellGlobalIndex = toRecvValuesByProc[i][j];
				//Store values
				RequestedValues[cellGlobalIndex] = std::vector<double>(nVariables, 0);
				for (int k = 0; k < nVariables; k++) {
					RequestedValues[cellGlobalIndex][k] = recvbuf[pointer+k];
				};				
				//Increment pointer
				pointer += nVariables;
			};			
		};
		
	};
};

#endif