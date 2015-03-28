#ifndef TURBO_MAINSOLVER_PerfomanceHelper
#define TURBO_MAINSOLVER_PerfomanceHelper

#include <cassert>
#include <chrono>
#include <mutex>
#include <memory>
#include <thread>
#include "parallelHelper.h"
#include "Timer.h"

//Structure representing hierarchy of phases
class PerfomancePhase {
	std::shared_ptr<ParallelHelper> _parallelHelper;
	std::mutex _isActiveMutex;
	std::chrono::time_point<std::chrono::high_resolution_clock> _startTime;
	std::chrono::time_point<std::chrono::high_resolution_clock> _endTime;
	std::chrono::time_point<std::chrono::high_resolution_clock> _beforeSyncTime;
	std::chrono::time_point<std::chrono::high_resolution_clock> _afterSyncTime;
	std::chrono::high_resolution_clock::duration _idleDuration;
	std::chrono::high_resolution_clock::duration _activeDuration;
	std::chrono::high_resolution_clock::duration _totalDuration;

	//Make class non copyable TO DO c++11 style
	//PerfomancePhase(const PerfomancePhase&) = delete;
    //PerfomancePhase& operator=(const PerfomancePhase&) = delete;
	PerfomancePhase(const PerfomancePhase&) {}; // non-copyable ('= delete' in C++11)
	PerfomancePhase& operator=(const PerfomancePhase&) {}; // non-copyable ('= delete' in C++11)

	//Time conversion functions
	template <class PeriodType> 
	inline double _getIdleTime() {
		//using durationType = std::chrono::duration<double, PeriodType>;
		//return std::chrono::duration_cast<durationType>(_idleDuration).count();
		return std::chrono::duration_cast<std::chrono::duration<double, PeriodType> >(_idleDuration).count();
	};

	template <class PeriodType> 
	inline double _getActiveTime() {
		//using durationType = std::chrono::duration<double, PeriodType>;
		//return std::chrono::duration_cast<durationType>(_activeDuration).count();
		return std::chrono::duration_cast<std::chrono::duration<double, PeriodType> >(_activeDuration).count();
	};

public:
	//Public fields
	std::string PhaseName;
	double IdleTime; //ms
	double ActiveTime; //ms
	double IsActive;
	std::map<std::string, std::shared_ptr<PerfomancePhase>> subPhases;

	//Constructor
	PerfomancePhase(std::shared_ptr<ParallelHelper>& parallelHelper, std::string phaseName) {
		PhaseName = phaseName;
		IsActive = false;
		_parallelHelper = parallelHelper;
		subPhases.clear();
	};

	//Add subphase
	std::weak_ptr<PerfomancePhase> CreateSubPhase(std::string subPhaseName) {
		std::shared_ptr<PerfomancePhase> newPhase(new PerfomancePhase(_parallelHelper, subPhaseName));
		subPhases[subPhaseName] = std::move(newPhase);
		return subPhases[subPhaseName];
	};
	
	//Start phase timer
	void Start() {
		if (_isActiveMutex.try_lock()) {
			_startTime = std::chrono::high_resolution_clock::now();
			IsActive = true;
			_idleDuration = _parallelHelper->getIdleTime();
			_totalDuration = std::chrono::high_resolution_clock::duration(0);
			subPhases.clear();
			_isActiveMutex.unlock();
		};
	};

	//Stop phase timer
	void Stop() {
		//Stop every unfinished subphase
		for (auto& subPhase : subPhases) subPhase.second->Stop(); 

		//Wait for all the threads and stop
		if (IsActive) {
			_isActiveMutex.lock(); // multithread sync barrier point
			_endTime = std::chrono::high_resolution_clock::now();
			IsActive = false;
			_isActiveMutex.unlock();
		};

		// compute idle duration based on parallel helper measurements
		_idleDuration = _parallelHelper->getIdleTime() - _idleDuration;

		// compute active stage duration
		_totalDuration = _endTime - _startTime;
		_activeDuration = _totalDuration - _idleDuration;
	};

	//Sync across MPI nodes
	void Sync() {
		_parallelHelper->Barrier(); //MPI sync barrier point
	};

	//Stop phase timer and sync across MPI nodes
	void StopAndSync() {
		Sync(); // sync
		Stop(); // stop
	};

	//Accessing perfomance information
	inline double GetIdleTimeMilliseconds() {
		double time = 0;
		for (auto& kvPair : subPhases) time += kvPair.second->_getIdleTime<std::milli>();
		time += _getIdleTime<std::milli>();
		return time;
	};

	inline double GetActiveTimeMilliseconds() {
		double time = _getActiveTime<std::milli>();
		return time;
	};

	inline double GetTotalTimeMilliseconds() {
		double time = _getActiveTime<std::milli>();
		return time;
	};
};

//Class that implements perfomance related functions
class PerfomanceHelper {
	std::shared_ptr<ParallelHelper> _parallelHelper; // reference to parallel part
public:
	std::unique_ptr<PerfomancePhase> RootPhase; //root phase that represents whole program
	std::shared_ptr<PerfomancePhase> TimeStepPhase; // reference to time step phase (must be set)
	Timer Timer; //timer for convenience usage

	//Constructor
	PerfomanceHelper() {
	};

	//Initialize
	void Init(ParallelHelper& parallelHelper) {
		//Initialize reference to parallel helper
		_parallelHelper = std::shared_ptr<ParallelHelper>(&parallelHelper);

		//Create and start root phase
		RootPhase = std::move(std::unique_ptr<PerfomancePhase>(new PerfomancePhase(_parallelHelper, "RootPhase")));
		RootPhase->Start();
	};

	//Finalize
	void Finalize() {
		//Stop and sync root phase
		RootPhase->StopAndSync();
	};

	//Initialize phases structure
	void InitPhases() {
	};

	//Save perfomance history
	void SavePerfomanceHistory() {

	};

	//Usefull functions that implement timer functionality
	std::chrono::high_resolution_clock::duration GetElapsedTime() {
		auto elapsedTime = Timer.ElapsedTime;
		return elapsedTime;
	};
};

#endif