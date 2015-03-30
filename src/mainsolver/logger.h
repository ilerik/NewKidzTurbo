#ifndef TURBO_MAINSOLVER_LOGGER
#define TURBO_MAINSOLVER_LOGGER

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include "ParallelManager.h"

//All possible logger mesage types
enum class LoggerMessageType {
	Error = 0,
	Warning = 1,
	Information = 2
};

//Possible message levels
enum class LoggerMessageLevel {
	Global = 0,
	Local = 1
};

//Implementation of logging mechanism
class Logger {
private:
	Logger& operator=(const Logger& logger); //non copyable

	//MPI parallel implementation
	std::shared_ptr<ParallelManager> _MPIManager;

	//Message level
	LoggerMessageLevel _level;

	//Output stream
	std::shared_ptr<std::ostream> _outputStream;
public:		
	//Constructor opens file stream
	Logger(std::shared_ptr<ParallelManager>& MPIManager, LoggerMessageLevel level, std::string filename) : _MPIManager(MPIManager), _level(level) {
		//Global level logger outputs only master node messages
		if ((_level == LoggerMessageLevel::Global) && (_MPIManager->IsMaster())) {
			std::ostringstream fname;
			fname<<filename<<".log";
			_outputStream = std::shared_ptr<std::ostream>(new std::ofstream(fname.str(), std::ios_base::out));
		};
		//Local level logger outputs node messages
		if ((_level == LoggerMessageLevel::Local)) {
			std::ostringstream fname;
			fname<<filename<<_MPIManager->rank()<<".log";
			_outputStream = std::shared_ptr<std::ostream>(new std::ofstream(fname.str(), std::ios_base::out));
		};
	};

	//Constructor uses specified stream
	Logger(std::shared_ptr<ParallelManager>& MPIManager, LoggerMessageLevel level, std::shared_ptr<std::ostream> outputStream) : _MPIManager(MPIManager), _level(level), _outputStream(outputStream) { 	};

	//Destructor flushes stream
	~Logger() {
		_outputStream->flush();
	};

	//Output stream style 
	template <typename T>
	Logger& operator<<(T const& value) {
		//Global level logger outputs only master node messages
		if ((_level == LoggerMessageLevel::Global) && (_MPIManager->IsMaster())) {
			logger._outputStream << value;
			return log;
		};
		logger._outputStream << value;
		return log;
	};
	
};

#endif