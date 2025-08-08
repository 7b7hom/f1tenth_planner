#include "graph_planner.hpp"

unique_ptr<string> Load(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open INI file: " << filename << endl;
    }

    string line;
    bool in_section = false;
    while (getline(file, line)) {
        // line이라는 문자열에서 해당 문자열을 찾지 못하면 npos를 반환
        if (line.find("[DRIVING_TASK]") != string::npos) {
            in_section = true;
            continue;
        }

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
    throw logic_error("Unreachable exit of the while loop!");
}

// Debug용 함수: map의 columns, rows 개수 print  
void map_size(DMap& map) {
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label);
}

void addDVectorToMap(DMap &map, string attr) {

    int len = 0;
    if (!map.empty()) {
        len = static_cast<int>(map.begin()->second.size());
    } else {
        throw invalid_argument("Empty Map! - addDVectorToMap");
    }

    DVector x_out(len), y_out(len);
    string x_label = "x_" + attr;
    string y_label = "y_" + attr;
    
    if (!attr.compare("bound_r")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[POS_X][i] + map[NORM_X][i] * map[WIDTH_R][i];
            y_out[i] = map[POS_Y][i] + map[NORM_Y][i] * map[WIDTH_R][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("bound_l")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[POS_X][i] - map[NORM_X][i] * map[WIDTH_L][i];
            y_out[i] = map[POS_Y][i] - map[NORM_Y][i] * map[WIDTH_L][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("raceline")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[POS_X][i] + map[NORM_X][i] * map[NORM_L][i];
            y_out[i] = map[POS_Y][i] + map[NORM_Y][i] * map[NORM_L][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    // i번째와 i-1번째 point의 delta_s 계산 
    // delta_s[0] = 0 
    else if (!attr.compare("delta_s")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len - 1; ++i) {
            x_out[i] = map[RL_S][i+1] - map[RL_S][i]; // 마지막 원소는 0
        }
        map[attr] = x_out; 
    }

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

void calcHeading(DVector &x_raceline,
                 DVector &y_raceline,
                 DVector &psi) {

    size_t N = x_raceline.size();
    psi.resize(N);

    // 닫힌 회로 가정. 예외 처리 필요
    double dx, dy;
    for (size_t i = 0; i < N; ++i) {
        
        if (i != N -1) {
            dx = x_raceline[i+1] - x_raceline[i];
            dy = y_raceline[i+1] - y_raceline[i];
        } else {
            dx = x_raceline[0] - x_raceline[N - 1];
            dy = y_raceline[0] - y_raceline[N - 1];
        } 
    psi[i] = atan2(dy, dx) - M_PI_2;
        
    normalizeAngle(psi[i]);

    }
    // cout << i<< ": " << psi[i] << endl;
    // cout << psi.size() << endl;

}

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

bool checkInsideBounds(const Vector2d& pos, const float veh_width) {

    if (stMap.find(LB_X) == stMap.end() || 
    stMap.find(LB_Y) == stMap.end() ||
    stMap.find(RB_X) == stMap.end() || 
    stMap.find(RB_Y) == stMap.end()) {
    throw invalid_argument("Boundary keys are missing in stMap!");
}

    int n = stMap[LB_X].size();
    MatrixXd bound_l(n,2);
    MatrixXd bound_r(n,2);
    for (int i = 0; i < n; ++i) {
        bound_l(i, 0) = stMap[LB_X][i];
        bound_l(i, 1) = stMap[LB_Y][i];

        bound_r(i, 0) = stMap[RB_X][i];
        bound_r(i, 1) = stMap[RB_Y][i];
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


    // cout << "-------here" << endl;
    // cout << dist_to_left_bound << endl;
    // cout << dist_to_right_bound << endl;
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