#ifndef TURBO_MAINSOLVER_PARALLELHELPER
#define TURBO_MAINSOLVER_PARALLELHELPER

#include <mpi.h>
#include <cassert>
#include <chrono>
#include "GasModel.h"

//Class that implements MPI-based operations for distributed memory architecture
class ParallelHelper {
	std::chrono::high_resolution_clock::duration _idleDuration; //total idle time
	bool isInitilized; //was MPI initialization successful
	bool _isMain; //was MPI initialization successful
	int _nProcessors; //total number of processes
	int _rank;  //current processor rank in _comm
	MPI_Comm _comm;	//global communicator
	std::vector<GasModel* > _gasModels; //References to gas models in use
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

	inline std::chrono::high_resolution_clock::duration getIdleTime() {
		return _idleDuration;
	};

	//Constructor
	ParallelHelper() {
		_isMain = false;
		isInitilized = false;
	};

	//Initialize MPI programm
	void Init(int *argc, char **argv[]) {
		MPI_Init(argc, argv);		
		_comm = MPI_COMM_WORLD;
		MPI_Comm_size(_comm, &_nProcessors);
		MPI_Comm_rank(_comm, &_rank);
		isInitilized = true;
		_isMain = true;
		_idleDuration = std::chrono::high_resolution_clock::duration(0);
	};

	void Init(MPI_Comm comm) {
		_comm = comm;
		int result = MPI_Comm_size(_comm, &_nProcessors);
		if (result != MPI_SUCCESS) throw new Exception("Init failed");
		MPI_Comm_rank(_comm, &_rank);
		isInitilized = true;
		_isMain = false;
		_idleDuration = std::chrono::high_resolution_clock::duration(0);
	};

	//Finilize MPI programm
	void Finalize() {
		if (_isMain) MPI_Finalize();
		isInitilized = false;
	};

	//Is master node
	inline bool IsMaster() {
		return _rank == 0;
	};

	//Barrier
	void Barrier() {
		auto before = std::chrono::high_resolution_clock::now();
		MPI_Barrier(_comm);
		auto after = std::chrono::high_resolution_clock::now();
		_idleDuration += after - before; // update idle time duration
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

	//Gather arrays of int of different size on master node
	void GathervInt(std::vector<cgsize_t>& local, std::vector<int>& counts, std::vector<cgsize_t>& result) {		
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
		if (result.size() != totalSize) result.resize(totalSize);
		MPI_Allgatherv(&local.front(), local.size(), MPI_INT, &result[0], &counts[0], &displs[0], MPI_INT, _comm);		
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

	//Maximum
	inline double Max(double x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_MAX, _comm);
		return res;
	};


public:
	//Problem oriented routines and data
	int nVariables; //Number of conservative variables
	std::vector<int> part; //Partitioning
	std::map<int, int> RequestedGasModelIndex; //Requested material index for cells 
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
	void InitExchange(std::vector<GasModel*>& gasModels) {
		//Initialize gas models
		_gasModels = gasModels;

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

	//Exchange gas model indexes of cells
	void ExchangeGasModelIndexes(Grid& grid, std::vector<int>& gasmodelIndexes) {

	};

	//Exchange values
	void ExchangeValues(Grid& grid, std::vector<double>& values) {
		//printf("00");
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
		//printf("11");

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