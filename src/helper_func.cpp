#include "graph_planner.hpp"

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

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

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


void addDVectorToMap(DMap &map,
                     string attr,
                     const IVector *idx_array) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();
    } 
    else {
        len = idx_array->size();
    }
    // cout << "attr: "<< attr << " / len:" << len << endl;

    DVector x_out(len), y_out(len);
    string x_label = "x_" + attr;
    string y_label = "y_" + attr;
    
    if (!attr.compare("bound_r")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__width_right][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__width_right][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("bound_l")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] - map[__x_normvec][i] * map[__width_left][i];
            y_out[i] = map[__y_ref][i] - map[__y_normvec][i] * map[__width_left][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("raceline")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__alpha][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__alpha][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    // i번째와 i-1번째 point의 delta_s 계산 
    // delta_s[0] = 0 
    else if (!attr.compare("delta_s")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len - 1; ++i) {
            x_out[i] = map[__s_racetraj][i+1] - map[__s_racetraj][i]; // 마지막 원소는 0
        }
        map[attr] = x_out; 
    }

    // map_size(map);
}

bool checkInsideBounds(const Vector2d& pos, const float veh_width) {

    if (sampling_map.find(__x_bound_l) == sampling_map.end() || 
    sampling_map.find(__y_bound_l) == sampling_map.end() ||
    sampling_map.find(__x_bound_r) == sampling_map.end() || 
    sampling_map.find(__y_bound_r) == sampling_map.end()) {
    throw invalid_argument("Boundary keys are missing in sampling_map!");
}

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

    double dist_to_left_bound = sqrt(d_bl_2);
    double dist_to_right_bound = sqrt(d_br_2);


    cout << "-------here" << endl;
    cout << dist_to_left_bound << endl;
    cout << dist_to_right_bound << endl;
    // VEH_WIDTH 조건 확인
    if (dist_to_left_bound < veh_width || dist_to_right_bound < veh_width)
    {
        // throw invalid_argument("Spline point violates VEH_WIDTH constraints!");
        return false;
    }

    // bound 밖에 있는지 여부 확인
    bool within_bounds = !(d_bl_2 > d_track2 || d_br_2 > d_track2);
    return within_bounds;
}

void printSplineInfo(const SplineMap& splineMap, const NodeMap& nodesPerLayer) {

    for (auto& [startPoint, endPoints] : splineMap) {
        for (auto& [endPoint, spline] : endPoints) {


        const Node& startNode = nodesPerLayer[startPoint.first][startPoint.second];
        const Node& endNode = nodesPerLayer[endPoint.first][endPoint.second];

        cout << "\n(" << startPoint.first << ", " << startPoint.second << ") --> ("
                  << endPoint.first << ", " << endPoint.second << ")\n";

        cout << "  [Start Node] x: " << startNode.x
                  << ", y: " << startNode.y
                  << ", psi: " << startNode.psi << "\n";
        cout << "  [End Node]   x: " << endNode.x
                  << ", y: " << endNode.y
                  << ", psi: " << endNode.psi << "\n";

        cout << "  coeffs_x (" << spline.coeffs_x.rows() << "x" << spline.coeffs_x.cols() << "):\n";
        cout << spline.coeffs_x << "\n";

        cout << "  coeffs_y (" << spline.coeffs_y.rows() << "x" << spline.coeffs_y.cols() << "):\n";
        cout << spline.coeffs_y << "\n";

        cout << "----------------------------------------";
        }
    }
}