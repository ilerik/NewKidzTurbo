#ifndef TURBO_MAINSOLVER_LOGGER
#define TURBO_MAINSOLVER_LOGGER

#include <string>

//Implementation of logging mechanism


class Logger {
public:
	int InitLogging(std::string filename) {
		return 0;
	};

	int WriteMessage(int priority, std::string msg) {
		
	};
};

#endif