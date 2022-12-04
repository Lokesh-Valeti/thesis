#ifndef __PREPROC_HPP__
#define __PREPROC_HPP__

#include "mpcio.hpp"

void preprocessing_comp(MPCIO &mpcio, int num_threads, char **args);
void preprocessing_server(MPCServerIO &mpcio, int num_threads, char **args);

#endif
