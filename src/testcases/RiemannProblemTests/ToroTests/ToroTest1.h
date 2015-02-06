#ifndef TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1
#define TURBO_TestCases_RiemannProblemTests_ToroTests_ToroTest1

#include "ToroTest.h"

namespace ToroTests {

//classes of different solvers for Test 1 (E.Toro's book // parag. 4.3.3, pp. 129 - 133)
class ToroTest1 : public ToroTest {
public:
	
	//constructor
	ToroTest1() {
		_riemannSolverType = RiemannSolverConfiguration::RiemannSolverType::Roe;
	};
}; //ToroTest1

}; //namespace


#endif