#include <vector>

#ifndef INCLUDE_HEATBATH_HH_
#define INCLUDE_HEATBATH_HH_

void do_sweep(double* gauge_field, int T, int L, double beta, const std::vector<int>& bound_ts);

#endif /* INCLUDE_HEATBATH_HH_ */
