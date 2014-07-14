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
			MPI_Gatherv(&local[0], local.size(), MPI_INT, &result[0], &counts[0], &displs[0], MPI_INT, 0, _comm);
		} else {
			MPI_Gatherv(&local[0], local.size(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, _comm);
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
			MPI_Gatherv(&local[0], local.size(), MPI_DOUBLE, &result[0], &counts[0], &displs[0], MPI_DOUBLE, 0, _comm);
		} else {
			MPI_Gatherv(&local[0], local.size(), MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, _comm);
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
};

#endif