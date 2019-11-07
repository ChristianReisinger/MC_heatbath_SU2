// io.hh

#ifndef __IO_HH__

#define __IO_HH__

// ********************

// Reads/writes a gauge field configuration.

// for all it ...
//   for all ix ...
//     for all iy ...
//       for all iz ...
//         for all mu ...

// Data starts after "BEGIN".

// h[0], h[1], h[2], h[3] for each link.

void read_gauge_field(double* const gauge_field, const char* filename, const int T, const int L);

void write_gauge_field(const double* const gauge_field, const char* filename, const int T, const int L, const char* header);

#endif
