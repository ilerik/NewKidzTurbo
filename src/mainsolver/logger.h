#ifndef TURBO_MAINSOLVER_LOGGER
#define TURBO_MAINSOLVER_LOGGER

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "parallelHelper.h"

//All possible logger mesage types
enum class LoggerMessageType {
	FATAL_ERROR = 0,
	INFORMATION = 1
};

//Possible message levels
enum class LoggerMessageLevel {
	GLOBAL = 0,
	LOCAL = 1
};

//Implementation of logging mechanism
class Logger {
private:
	Logger& operator=(const Logger& logger); //non copyable
	int _nProcessors;
	int _rank;
	ParallelHelper* _parallelHelper;
	std::ofstream _output;
public:		
	int InitLogging(ParallelHelper& parallelHelper, std::string filename) {
		_parallelHelper = &parallelHelper;
		_nProcessors = _parallelHelper->getProcessorNumber();
		_rank = _parallelHelper->getRank();		
		std::ostringstream fname;
		fname<<filename<<_rank;
		_output.open(fname.str(), std::ios_base::out);
		return 0;
	};

	int WriteMessage(LoggerMessageLevel level, LoggerMessageType type, std::string msg) {
		if (level == LoggerMessageLevel::GLOBAL) {
			if (_rank == 0) {
				_output<<"GLOBAL : "<<msg<<"\n";
				_output.flush();
				//std::cout<<"GLOBAL : "<<msg<<"\n";
			};			
		};
		if (level == LoggerMessageLevel::LOCAL) {						
			_output<<"LOCAL rank = "<<_rank<<" : "<<msg<<"\n";
			_output.flush();			
			//std::cout<<"LOCAL  "<<_rank<<" : "<<msg<<"\n";
		};
		return 0;
	};

	int WriteMessage(LoggerMessageLevel level, LoggerMessageType type, std::string msg, int n) {
		if (level == LoggerMessageLevel::GLOBAL) {
			if (_rank == 0) {
				_output<<"GLOBAL : "<<msg<<n<<"\n";
				_output.flush();				
			};			
		};
		if (level == LoggerMessageLevel::LOCAL) {						
			_output<<"LOCAL rank = "<<_rank<<" : "<<msg<<n<<"\n";
			_output.flush();						
		};
		return 0;
	};

	int FinilizeLogging() {
		_output.close();
		return 0;
	};
	
};

#endif