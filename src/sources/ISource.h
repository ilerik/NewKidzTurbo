#ifndef NewKidzTurbo_Sources_ISource
#define NewKidzTurbo_Sources_ISource

#include "cgnslib.h"
#include "configuration.h"

//Source term base class
class ISource {
public:	
	virtual std::vector<double> GetResidual(std::vector<double>& values) = 0;
};

#endif
