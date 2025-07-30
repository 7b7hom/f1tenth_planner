#include "graph_planner.hpp"

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label);
}

// DMap을 CSV에 작성 
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter) {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();

    // Header
    size_t i = 0;
    for (const auto& [key, _] : map) {
        file << key;
        if (++i != num_cols) file << delimiter;
    }
    file << '\n';

    // Row map
    for (size_t row = 0; row < num_rows; ++row) {
        size_t j = 0;
        for (const auto& [_, col] : map) {
            file << col[row];
            if (++j != num_cols) file << delimiter;
        }
        file << '\n';
    }

    file.close();
}

// Debug용 함수: map의 columns, rows 개수 print  
void map_size(DMap& map) {
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

unique_ptr<string> Load(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open INI file: " << filename << endl;
    }

    string line;
    bool in_section = false;
    while (getline(file, line)) {
        // 섹션 시작
        if (line.find("[DRIVING_TASK]") != string::npos) {
            in_section = true;
            continue;
        }

        // 다른 섹션으로 넘어가면 종료
        if (in_section && line.find('[') != string::npos)
            break;

        // track 키 찾기
        if (in_section && line.find("track") != string::npos) {
            size_t eq_pos = line.find('=');
            if (eq_pos != string::npos) {
                string value = line.substr(eq_pos + 1);
                value.erase(0, value.find_first_not_of(" \t\r\n"));
                value.erase(value.find_last_not_of(" \t\r\n") + 1);
                return make_unique<string>(value);
            }
        }
    }

}

bool checkInsideBounds(const Vector2d& pos) {
    int n = sampling_map[__x_bound_l].size();
    MatrixXd bound_l(n,2);
    MatrixXd bound_r(n,2);
    for (int i = 0; i < n; ++i) {
        bound_l(i, 0) = sampling_map[__x_bound_l][i];
        bound_l(i, 1) = sampling_map[__y_bound_l][i];

        bound_r(i, 0) = sampling_map[__x_bound_r][i];
        bound_r(i, 1) = sampling_map[__y_bound_r][i];
    }
    
    MatrixXd centerline = (bound_l + bound_r) / 2;

    // 가장 가까운 segment 인덱스 찾기
    int closest_idx = -1;
    double min_dist2 = numeric_limits<double>::max();
    for (int i = 0; i < centerline.rows() - 1; ++i) {
        // segment 중심 계산
        Vector2d mid = (centerline.row(i) + centerline.row(i + 1)) / 2.0;
        double dist2 = (mid - pos).squaredNorm();
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closest_idx = i;
        }
    }

    if (closest_idx < 0 || closest_idx >= bound_l.rows() - 1)
        return false; // 예외 처리

    // bound_l, bound_r, centerline 보간 (선형 보간 10개 지점)
    int interp_points = 10;
    MatrixXd bl_interp(interp_points, 2);
    MatrixXd br_interp(interp_points, 2);
    MatrixXd center_interp(interp_points, 2);

    for (int i = 0; i < interp_points; ++i) {
        double t = static_cast<double>(i) / (interp_points - 1);
        bl_interp.row(i) = (1 - t) * bound_l.row(closest_idx) + t * bound_l.row(closest_idx + 1);
        br_interp.row(i) = (1 - t) * bound_r.row(closest_idx) + t * bound_r.row(closest_idx + 1);
        center_interp.row(i) = (1 - t) * centerline.row(closest_idx) + t * centerline.row(closest_idx + 1);
    }

    // pos에 가장 가까운 center_interp 인덱스 찾기
    int nearest_idx = -1;
    double best_dist2 = numeric_limits<double>::max();
    for (int i = 0; i < interp_points; ++i) {
        double d2 = (center_interp.row(i) - pos.transpose()).squaredNorm();
        if (d2 < best_dist2) {
            best_dist2 = d2;
            nearest_idx = i;
        }
    }

    // bound 사이 거리 (제곱)
    double d_track2 = (bl_interp.row(nearest_idx) - br_interp.row(nearest_idx)).squaredNorm();

    // 차량에서 각 bound까지 거리 (제곱)
    double d_bl_2 = (bl_interp.row(nearest_idx) - pos.transpose()).squaredNorm();
    double d_br_2 = (br_interp.row(nearest_idx) - pos.transpose()).squaredNorm();

    // bound 밖에 있는지 여부 확인
    bool within_bounds = !(d_bl_2 > d_track2 || d_br_2 > d_track2);
    return within_bounds;
}
