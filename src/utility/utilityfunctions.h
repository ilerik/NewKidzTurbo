#ifndef TURBO_UtilityFunctions
#define TURBO_UtilityFunctions

//Signum function from http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif