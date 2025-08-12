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
DMap readDMapFromCSV(const string& pathname) {
    DMap gtMap;
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels) 
        gtMap[label] = csv.GetColumn<double>(label);

    return gtMap;
}

pair<DVector, DVector> computeBoundRight(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &width_r) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || width_r.empty()) {
        throw runtime_error("computeBoundRight() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_bound_r(len), y_bound_r(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_bound_r[i] = pos_x[i] + norm_x[i] * width_r[i];
        y_bound_r[i] = pos_y[i] + norm_y[i] * width_r[i];
    }
    
    return {x_bound_r, y_bound_r};

}

pair<DVector, DVector> computeBoundLeft(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &width_l) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || width_l.empty()) {
        throw runtime_error("computeBoundLeft() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_bound_l(len), y_bound_l(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_bound_l[i] = pos_x[i] - norm_x[i] * width_l[i];
        y_bound_l[i] = pos_y[i] - norm_y[i] * width_l[i];
    }
    
    return {x_bound_l, y_bound_l};

}

pair<DVector, DVector> computeRaceline(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &norm_l) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || norm_l.empty()) {
        throw runtime_error("computeBoundRaceline() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_raceline(len), y_raceline(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_raceline[i] = pos_x[i] + norm_x[i] * norm_l[i];
        y_raceline[i] = pos_y[i] + norm_y[i] * norm_l[i];
    }
    
    return {x_raceline, y_raceline};

}

DVector computeDeltaS(DVector &rl_s) {
    if (rl_s.empty()) {
        throw runtime_error("computeDeltaS() - Empty DVector !!");
    }

    int len = rl_s.size();
    DVector rl_ds(len);

    // 마지막 원소는 0
    for (size_t i = 0; i < len - 1; ++i) {
        rl_ds[i] = rl_s[i+1] - rl_s[i];
    }
    
    return rl_ds;

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

DVector calcHeading(DVector &x_raceline, DVector &y_raceline) {

    DVector psi;
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
    return psi;
}

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
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