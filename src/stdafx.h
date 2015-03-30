//Boost library files
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/algorithm/string/replace.hpp>

//Usefull macro

//CGNS error handling macro
#define CALL_CGNS(cgns_func) {                                                 \
                               const int ierr = cgns_func;                     \
                               if (ierr)                                       \
                               {                                               \
                                 const char * error_msg = cg_get_error();      \
                                 cg_error_exit();                              \
                               }                                               \
                             };