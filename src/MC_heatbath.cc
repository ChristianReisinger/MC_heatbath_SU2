// MC_heatbath.cc

// I do not remember, where I got these numbers --> be cautious!
// beta = 2.3   -->   a = 0.171 fm
// beta = 2.4   -->   a = 0.128 fm
// beta = 2.5   -->   a = 0.089 fm

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <getopt.h>
#include <vector>
#include <set>

#include <fields.hh>
#include <geometry.hh>
#include <io.hh>
#include <linear_algebra.hh>
#include <ranlux.hh>
#include <Wilson_loops.hh>
#include <heatbath.hh>

#include <helper_functions.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace MC_heatbath_SU2 {

void handle_GNU_options(int argc, char **&argv, std::set<int> &bound_ts, std::string &prefix) {

	static struct option long_opts[] = {
			{ "fixed-bounds", required_argument, 0, 'f' },
			{ "prefix", required_argument, 0, 'p' },
			{ 0, 0, 0, 0 }
	};
	int opt = -1, long_opts_i = 0;
	while ((opt = getopt_long(argc, argv, "f:p:", long_opts, &long_opts_i)) != -1) {
		switch (opt) {
			case 'f': {
				std::vector<int> bound_ts_vec = parse_unsigned_int_list(optarg);
				bound_ts = std::set<int>(bound_ts_vec.begin(), bound_ts_vec.end());
			}
			break;
			case 'p':
				prefix = optarg;
			break;
		}
	}
	argv = argv + optind - 1;

}

}
}
}

int main(int argc, char **argv) {
	using namespace std;
	using namespace de_uni_frankfurt_itp::reisinger::MC_heatbath_SU2;

	if (argc < 9 || argc > 14) {
		cout << "Usage: " << argv[0]
				<< " [(--fixed-bounds | -f) <ts>] [(--prefix | -p) <prfx>] <output_dir> <beta> <T> <L> <seed> (cold | hot | <init_file>) <num_MC_sweeps_max> <num_MC_sweeps_out>\n";
		return 0;
	}

	set<int> bound_ts;
	string prefix("conf");
	handle_GNU_options(argc, argv, bound_ts, prefix);

	string path(argv[1]);
	double beta;
	int T, L, seed, num_MC_sweeps_max, num_MC_sweeps_out;

	stringstream args_ss;
	args_ss << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[7] << " " << argv[8];
	args_ss >> beta >> T >> L >> seed >> num_MC_sweeps_max >> num_MC_sweeps_out;

	if (beta <= 0.0) {
		cerr << "Error: beta <= 0.0\n";
		return 0;
	}

	if (T < 2 || L < 2) {
		cerr << "Error: T < 2 or L < 2\n";
		return 0;
	}

	if (seed < 1) {
		cerr << "Error: seed < 1\n";
		return 0;
	}

	if (num_MC_sweeps_max < 1 || num_MC_sweeps_out < 1) {
		cerr << "Error: num_MC_sweeps_max < 1 or num_MC_sweeps_out < 1.\n";
		return 0;
	}

	int config_id_digits = floor(log10(num_MC_sweeps_max)) + 1;

	InitializeRand(seed);

	double *gauge_field;
	Gauge_Field_Alloc(&gauge_field, T, L);

	if (string(argv[6]) == string("hot")) {
		cout << "Initializing hot start ...\n";
		Gauge_Field_Random(gauge_field, T, L);
	} else if (string(argv[6]) == string("cold")) {
		cout << "Initializing cold start ...\n";
		Gauge_Field_Unity(gauge_field, T, L);
	} else {
		cout << "Initializing config from file ...\n";
		read_gauge_field(gauge_field, argv[6], T, L);
	}

	double P = Average_Plaquette(gauge_field, T, L);
	cout << "  ... <P> = " << P << ".\n";

	// Update links using the heatbath algorithm.

	ostringstream plaq_filename_oss;
	plaq_filename_oss << path << "/P.dat";
	ofstream plaquette_ofs(plaq_filename_oss.str());
	if (plaquette_ofs.fail()) {
		cerr << "Error: Unable to open average plaquette file.\n";
		return 0;
	}

	for (int MC_sweep_i = 1; MC_sweep_i <= num_MC_sweeps_max; MC_sweep_i++) {
		do_sweep(gauge_field, T, L, beta, bound_ts);
		cout << "  ... sweep " << MC_sweep_i << " - end ...\n";

		P = Average_Plaquette(gauge_field, T, L);
		cout << "  ... <P> = " << P << ".\n";
		plaquette_ofs << setw(5) << MC_sweep_i << " " << showpos << setprecision(5) << scientific << P << noshowpos << endl;

		if (MC_sweep_i % num_MC_sweeps_out == 0) {
			ostringstream config_filename_oss;
			config_filename_oss << path << "/" << prefix << "." << setfill('0') << setw(config_id_digits) << MC_sweep_i;

			ostringstream header_oss(argv[0]);
			for (int arg_i = 1; arg_i < argc; ++arg_i)
				header_oss << " " << argv[arg_i];
			write_gauge_field(gauge_field, config_filename_oss.str().c_str(), T, L, header_oss.str().c_str());
		}
	}

	Gauge_Field_Free(&gauge_field);

	return EXIT_SUCCESS;
}

// ********************
