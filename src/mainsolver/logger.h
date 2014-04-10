#ifndef TURBO_MAINSOLVER_LOGGER
#define TURBO_MAINSOLVER_LOGGER

#include <string>
#include <iostream>
#include <sstream>

//All possible logger mesage types
enum LoggerMessageType {
	FATAL_ERROR = 0,
	INFORMATION = 1
};

//Possible message levels
enum LoggerMessageLevel {
	GLOBAL = 0,
	LOCAL = 1
};

//Implementation of logging mechanism
class Logger {
	int _rank;
public:	
	int InitLogging(std::string filename, int rank) {
		_rank = rank;
		return 0;
	};

	int WriteMessage(LoggerMessageLevel level, LoggerMessageType type, std::string msg) {
		if (level == GLOBAL) {
			if (_rank == 0) std::cout<<"GLOBAL : "<<msg<<"\n";
		};
		if (level == LOCAL) {
			std::cout<<"LOCAL  "<<_rank<<" : "<<msg<<"\n";
		};
		return 0;
	};
};

#endif