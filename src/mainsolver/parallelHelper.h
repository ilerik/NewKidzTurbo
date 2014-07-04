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

};

#endif