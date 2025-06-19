#include "config.h"

void samplePointsFromRaceline(const vector<double> &kappa,
                              const vector<double> &dist,
                              double d_curve,
                              double d_straight,
                              double curve_th,
                              vector<int> &idx_array) {

    double next_dist = 0.0;
    double next_dist_min = 0.0;
    double cur_dist = 0.0;

    for (size_t idx = 0; idx < kappa.size(); ++idx) {
        double kappa_val = kappa[idx];
        double dist_val = dist[idx];

        if ((cur_dist + dist_val) > next_dist_min) {
            if (fabs(kappa_val) > curve_th) {
                next_dist = cur_dist;
            }
        }

        if ((cur_dist + dist_val) > next_dist) {
            idx_array.push_back(static_cast<int>(idx));

            if (fabs(kappa_val) < curve_th) { // 직선 구간
                next_dist += d_straight;
            } else { // 곡선 구간
                next_dist += d_curve;
            }
            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist_val;
    }

    ofstream file("inputs/idx_array.txt");
    for (int idx : idx_array) 
        file << idx << "\n";
}

int main() {
    Offline_Params params;
    Document csv("inputs/gtpl_levine.csv",
                LabelParams(0, -1),
                SeparatorParams(';'));

	setTrackData(csv);
    vector<int> idx_array;

    samplePointsFromRaceline(multiColumns[__kappa], 
                            multiColumns[__delta_s], 
                            params.LON_CURVE_STEP, 
                            params.LON_STRAIGHT_STEP,
                            params.CURVE_THR,
                            idx_array);

 return 0;
}