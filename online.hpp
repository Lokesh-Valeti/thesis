#ifndef __ONLINE_HPP__
#define __ONLINE_HPP__

#include "mpcio.hpp"

void online_comp(MPCIO &mpcio, int num_threads, char **args);
void online_server(MPCServerIO &mpcio, int num_threads, char **args);

#endif
