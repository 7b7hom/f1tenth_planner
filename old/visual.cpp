#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "setTrackData.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

bool readXYcsv(const string& filename, vector<double>& xs, vector<double>& ys)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "파일을 열 수 없습니다: " << filename << endl;
        return false;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string x_str, y_str;
        if (!getline(ss, x_str, ',')) continue;
        if (!getline(ss, y_str, ',')) continue;

        xs.push_back(stod(x_str));
        ys.push_back(stod(y_str));
    }
    return true;
}

int main(int argc, char *argv[])
{
	Document csv("inputs/gtpl_levine.csv",
				LabelParams(0, -1),
				SeparatorParams(';'));

	setTrackData(csv);

	ifstream file("inputs/idx_array.txt");
	int idx;

	if (file.is_open()) {
		while (file >> idx) {
			multiColumns[__x_sampling].push_back(multiColumns[__x_raceline][idx]);
			multiColumns[__y_sampling].push_back(multiColumns[__y_raceline][idx]);
		}
		file.close();
	} else {
		cerr << "파일 열기 실패" << endl;
	}

	plt::clf();

	plt::plot(multiColumns[__x_bound_l], multiColumns[__y_bound_l], {{"color", "black"}});
	plt::plot(multiColumns[__x_bound_r], multiColumns[__y_bound_r], {{"color", "black"}});
    plt::plot(multiColumns[__x_ref], multiColumns[__y_ref], {{"color", "blue"}});
    plt::plot(multiColumns[__x_raceline], multiColumns[__y_raceline], {{"color", "red"}});
    plt::scatter(multiColumns[__x_sampling], multiColumns[__y_sampling], 40.0, {{"color", "pink"}});

	plt::title("Track");
    plt::grid(true);
	plt::axis("equal");
	plt::show();

	return 0;
}





#if 0
int main()
{
	int n = 1000;
	std::vector<double> x, y, z;

	for(int i=0; i<n; i++) {
		x.push_back(i*i);
		y.push_back(sin(2*M_PI*i/360.0));
		z.push_back(log(i));

		if (i % 10 == 0) {
			// Clear previous plot
			plt::clf();
			// Plot line from given x and y data. Color is selected automatically.
			plt::plot(x, y);
			// Plot a line whose name will show up as "log(x)" in the legend.
			plt::named_plot("log(x)", x, z);

			// Set x-axis to interval [0,1000000]
			plt::xlim(0, n*n);

			// Add graph title
			plt::title("Sample figure");
			// Enable legend.
			plt::legend();
			// Display plot continuously
			plt::pause(0.01);
		}
	}
}
#endif 