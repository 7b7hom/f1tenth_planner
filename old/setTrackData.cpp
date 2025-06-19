#include "setTrackData.hpp"

void computeOffsetXY(
    const string& alpha, int sign,
    const string& out_x_name, const string& out_y_name
) {
    size_t N = multiColumns.at(__x_ref).size();
    vector<double> x_out(N), y_out(N);

    for (size_t i = 0; i < N; ++i) {
        x_out[i] = multiColumns.at(__x_ref)[i] + sign * multiColumns.at(__x_normvec)[i] * multiColumns.at(alpha)[i];
        y_out[i] = multiColumns.at(__y_ref)[i] + sign * multiColumns.at(__y_normvec)[i] * multiColumns.at(alpha)[i];
    }

    multiColumns[out_x_name] = x_out;
    multiColumns[out_y_name] = y_out;
}

void setTrackData(Document& csv) {

    readCSV(csv, multiColumns);

    computeOffsetXY(__width_right, +1, "x_bound_r", "y_bound_r");
    computeOffsetXY(__width_left, -1, "x_bound_l", "y_bound_l");
    computeOffsetXY(__alpha, +1, "x_raceline", "y_raceline");

    vector<double>& s_racetraj = multiColumns[__s_racetraj];
    // 구간별 s 저장     
    for (size_t i = 0; i < multiColumns[__s_racetraj].size()-1; ++i) {
        // if (i==0) delta_s.push_back(0);
        multiColumns[__delta_s].push_back(s_racetraj[i+1] - s_racetraj[i]);
    }

}
