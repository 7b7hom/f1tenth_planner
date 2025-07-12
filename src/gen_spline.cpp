#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label);
}


// DMap을 CSV에 작성 
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',') {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    size_t num_cols = map.size();   // 열 개수 = map의 key 개수
    size_t num_rows = map.begin()->second.size();   // 행 개수 = 첫 번째 key의 값의 길이

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

// Dvector를 Map 구조로 추가(연산)
void addDVectorToMap(DMap &map,
                     string attr,
                     const IVector *idx_array = nullptr) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();  // 지정X -> 전체 데이터 대상
    } 
    else {
        len = idx_array->size();    // 지정한 인덱스 수만큼
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

void samplePointsFromRaceline(const DVector& kappa,
                              const DVector& dist,
                              double d_curve,
                              double d_straight,
                              double curve_th,
                              IVector& idx_array) {

    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    for (size_t i = 0; i < n; ++i) {

        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist;
        }

        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) {
            idx_array.push_back(static_cast<int>(i));

            if (fabs(kappa[i]) < curve_th) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist[i];
    }

    // for (size_t i=0; i < idx_array.size(); ++i) 
    //     cout << idx_array[i] << endl;
    // cout << "size: " << idx_array.size() << endl;
}

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
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

void genNode(NodeMap& nodesPerLayer,
            IVector& raceline_index_array,
            const double veh_width,
            float lat_resolution) {
    
    const size_t N = sampling_map[__alpha].size();
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        node.layer_idx = i; 
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);
        
        // cout << "layer 길이" << (sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width/2)<< endl;
        // cout << "layer 내에서 raceline index:" << raceline_index << endl;
        // cout << "-----" << endl;

        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);    // 기준선에서의 위치
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);  // 기준선에서 수직한 노멀 벡터 따라 노드 배치
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지
        int node_idx = 0;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;  // num_nodes : 좌우 총 가능한 노드 수
        nodesPerLayer[i].resize(num_nodes); 

        // cout << i << "번째 layer의 node 개수는 " << num_nodes << endl;
        // node별 loop 
        for (double alpha = start_alpha; alpha <= sampling_map[__width_right][i] - veh_width / 2 ; alpha+=lat_resolution) {
            // node_alphas.push_back(alpha);
            // node의 좌표 계산.
            node_pos = ref_xy + alpha * norm_vec;
            // node의 layer내의 인덱스 계산.
            node.node_idx = node_idx;
            node.x = node_pos.x();
            node.y = node_pos.y();
            node.psi = 0.0;
            node.kappa = 0.0;        
            node.raceline = (node_idx == raceline_index);

            double psi_interp;
            if (node_idx < raceline_index) {
                
                if (abs(sampling_map[__psi_bound_l][i] - sampling_map[__psi][i]) >= M_PI) 
                {   
                    double bl = sampling_map[__psi_bound_l][i] + 2 * M_PI * (sampling_map[__psi_bound_l][i] < 0);
                    double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
                    psi_interp = bl + (p - bl) * node_idx / raceline_index;
                                       
                }
                else {
                    psi_interp = sampling_map[__psi_bound_l][i] + (sampling_map[__psi][i] - sampling_map[__psi_bound_l][i]) * (node_idx+1) / raceline_index;
                }
                node.psi = normalizeAngle(psi_interp);
            }
            else if (node_idx == raceline_index) {
                psi_interp = sampling_map[__psi][i];
                node.psi = psi_interp;
            }
            else {
                int remain = num_nodes - raceline_index - 1;
                double t = static_cast<double>(node_idx - raceline_index) / max(remain, 1);  // 0 ~ 1
                psi_interp = sampling_map[__psi][i] + t * (sampling_map[__psi_bound_r][i] - sampling_map[__psi][i]);
                node.psi = normalizeAngle(psi_interp);
            }
            // cout << i << "번째 레이어의" <<node_idx << "번째 노드의 psi는" << node.psi << endl;
            
            nodesPerLayer[i][node_idx] = node;
            ++node_idx;
        }

        // cout << i << "번째 Layer의" << endl;
        // for (size_t i =0; i < node_pos.size(); ++i) {        
        //     cout << i << "번째 Node" << endl;
        //     cout << node_pos[i] << endl;
        // }
                // 각 node의 psi, kappa 계산하는 로직 추가

    }
}

void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale = 0.5) {
    double dx, dy;
    double theta, arrow_len;
    double angle;
    double x_arrow1, y_arrow1;
    double x_arrow2, y_arrow2;

    for (size_t i = 0; i < x.size(); ++i) {
        dx = scale * cos(psi[i] + M_PI_2);
        dy = scale * sin(psi[i] + M_PI_2);
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        plt::plot(x_line, y_line, {{"color", "green"}});

        #if 1 // 화살촉 그리기
        theta = atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);
        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "green"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "green"}});
        #endif
    }
}

// NodeMap에 저장된 모든 노드들을 보라색 점으로 플로팅하고, 각 노드의 헤딩을 화살표로 시각화
void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            node_x.push_back(node.x);
            node_y.push_back(node.y);
        }
    }
    plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}, {"label", "Nodes"}});
}



// -------------------------------------------
// gen_spline --------------------------------
// -------------------------------------------

// 경로 행렬의 각 점 사이의 유클리드 거리 계산 함수
VectorXd computeEuclideanDistances(const MatrixXd& path) {
    int N = path.rows() - 1;    // 구간 개수는 점 개수 - 1
    VectorXd dists(N);          // 거리 저장 벡터
    for (int i = 0; i < N; ++i) {
        dists(i) = (path.row(i + 1) - path.row(i)).norm();  // 두 점 사이 거리 계산 (유클리드 노름)
    }
    return dists;
}

// path : x,y 좌표 배열
// el_lengths 전달 안 하면 nullptr (점 사이 거리는 필수 x)
// psi_s, psi_e도 전달 안 하면 NAN (필수 x)
SplineResult calcSplines(const MatrixXd& path, const VectorXd* el_lengths = nullptr,
                         double psi_s = NAN, double psi_e = NAN,
                         bool use_dist_scaling = true) {

        // 1.폐곡선 여부 판단
        MatrixXd updated_path = path;
        
        // 시작 점과 끝 점 사이 거리가 매우 작으면 같다고 간주한다.
        bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6;

        if (closed) {
            // 폐곡선 유지 : 시작점 다시 한 번 넣기 (마지막 스플라인 조각을 만들기 위해서)
            updated_path.conservativeResize(path.rows() + 1, Eigen::NoChange);
            updated_path.row(updated_path.rows() - 1) = path.row(0);
        }

        int no_splines = updated_path.rows() - 1;   // 경로 구간 수 = 점 개수 - 1

        VectorXd ds;
        // 유클리드 거리 담고 있는 el_lengths (nullptr -> 유클리드 거리 계산, 아니면 el_lengths가 가리키는 VerctorXd 객체 거리 사용)
        if (el_lengths == nullptr) {
            ds = computeEuclideanDistances(updated_path);  // 기본 거리 계산 함수 호출
        } else {
            ds = *el_lengths;  // 외부에서 받은 길이 벡터 사용
        }

        // 유클리드 거리 출력
        cout << "ds(거리벡터):\n" << ds.transpose() << "\n";

        // 2. 행렬 M 구성
        // 각 스플라인은 x, y 각각에 대해 4개의 계수(a0 ~ a3) 가짐
        // 스플라인 0 -> 행 0 ~ 3, 스플라인 1 -> 행 4 ~ 7

        // 위치 조건(각 스플라인은 시작점과 끝점에서 위치가 정확히 맞아야 함) : 식 2*(N-1)개
        // 1차 도함수 연속(기울기 연속) : 식 (N-2)개
        // 2차 도함수 연속(곡률 연속) : 식 (N-2)개
        // 양 끝 점(경계 조건) : 식 2개

        // 스케일링 없이 진행. 1, 2차 도합수 조건에 직접 거리 d 사용하면 됨
        
        MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
        VectorXd b_x = VectorXd::Zero(no_splines * 4);
        VectorXd b_y = VectorXd::Zero(no_splines * 4);

        // 위치 조건
        for (int i = 0; i < no_splines; i++) {
            int row_s = i * 4;  // 스플라인 시작 행
            double d = ds(i);

            // 시작점 위치 조건
            M(row_s, row_s) = 1.0;  // a_0 = p_i
            b_x(row_s) = updated_path(i,0);
            b_y(row_s) = updated_path(i,1);

            // 끝점 위치 조건
            RowVector4d T(1.0, d, d*d, d*d*d);
            M.block(row_s + 1, row_s, 1, 4) = T; // 시작 행, 시작 열, 선택할 행, 선택할 열
            b_x(row_s + 1) = updated_path(i+1, 0);
            b_y(row_s + 1) = updated_path(i+1, 1);
        }

        // 기울기 조건
        // i = 0 ~ no_splines - 2까지만 적용 (마지막 스플라인 뒤에는 비교할 스플라인이 없기 때문)
        for (int i = 0; i < no_splines - 1; ++i) {
            int row = 2 * no_splines + i;
            double d = ds(i);

            M(row, 4 * i + 1) = 1.0;              // a_1^i
            M(row, 4 * i + 2) = 2.0 * d;          // 2a_2^i * d
            M(row, 4 * i + 3) = 3.0 * d * d;      // 3a_3^i * d^2
            M(row, 4 * (i + 1) + 1) = -1.0;       // -a_1^{i+1}

            b_x(row) = 0.0;
            b_y(row) = 0.0;
        }

        // 곡률 조건
        for (int i = 0; i < no_splines - 1; ++i) {
            int row = 3 * no_splines - 1 + i;
            double d = ds(i);

            M(row, 4 * i + 2) = 2.0;                // 2a_2^i
            M(row, 4 * i + 3) = 6.0 * d;            // 6a_3^i * d
            M(row, 4 * (i + 1) + 2) = -2.0;         // -2a_2^{i+1}

            b_x(row) = 0.0;
            b_y(row) = 0.0;
        }

        // 닫히지 않은 경로에서 시작점과 끝점 방향 조건
        if (!closed) {
            // 시작점 도함수 조건 (t=0)
            int start_row = 4 * no_splines - 2; // 마지막에서 두 번째 줄이 시작점 도함수 조건
            M.row(start_row).setZero();
            M(start_row, 1) = 1.0;  // a_1항에 1 설정 (S'(0) = a_1)

            double d_s = ds(0); // 첫 번째 구간 길이
            b_x(start_row) = cos(psi_s) * d_s;  // x방향 단위벡터 * d_s
            b_y(start_row) = sin(psi_s) * d_s;  // y방향 단위벡터 * d_s

            // 끝점 도함수 조건 (t=d)
            int end_row = 4 * no_splines - 1;
            M.row(end_row).setZero();
            
            double d_e = ds(no_splines - 1);    // 마지막 구간 길이
            int offset = 4 * (no_splines - 1);  // 마지막 스플라인 계수들의 시작 인덱스
            M(end_row, offset + 1) = 1.0;
            M(end_row, offset + 2) = 2.0 * d_e;
            M(end_row, offset + 3) = 3.0 * d_e * d_e;

            b_x(end_row) = cos(psi_e) * d_e;
            b_y(end_row) = sin(psi_e) * d_e;
        }

        // [DEBUG] M, b_x, b_y 출력
        cout << "[DEBUG] M 행렬:\n" << M << "\n";
        cout << "[DEBUG] b_x 벡터:\n" << b_x.transpose() << "\n";
        cout << "[DEBUG] b_y 벡터:\n" << b_y.transpose() << "\n";

        // 연립방정식 풀기
        VectorXd x_les = M.colPivHouseholderQr().solve(b_x);
        VectorXd y_les = M.colPivHouseholderQr().solve(b_y);

        MatrixXd coeffs_x = Map<MatrixXd>(x_les.data(), 4, no_splines).transpose();
        MatrixXd coeffs_y = Map<MatrixXd>(y_les.data(), 4, no_splines).transpose();


        // 법선 벡터 계산 (노멀 벡터)
        MatrixXd normvec(no_splines, 2);
        for (int i = 0; i < no_splines; ++i) {
            double dx = coeffs_x(i, 1);  // a₁: x방향 도함수 시작점
            double dy = coeffs_y(i, 1);  // a₁: y방향 도함수 시작점
            normvec(i, 0) = -dy;
            normvec(i, 1) = dx;
        }

       // 정규화된 법선 벡터 계산
        MatrixXd normvec_normalized(no_splines, 2);
        for (int i = 0; i < no_splines; ++i) {
            double norm = normvec.row(i).norm();
            if (norm > 1e-6)
                normvec_normalized.row(i) = normvec.row(i) / norm;
            else
                normvec_normalized.row(i) = Vector2d(0.0, 0.0); // norm이 너무 작을 때 0으로 대체
        }

        // 결과 구조체에 저장
        SplineResult result;
        result.coeffs_x = coeffs_x;
        result.coeffs_y = coeffs_y;
        result.M = M;
        result.normvec_normalized = normvec_normalized;


        return result;

}

// ------------------evaluate spline------------------

struct SplineEval {
    Vector2d pos;    // 위치 (x, y)
    double heading;  // 헤딩 (rad)
    double curvature; // 곡률 (1/m)
};

SplineEval evaluateSpline(const RowVector4d& coeffs_x,
                          const RowVector4d& coeffs_y,
                          double t) {
    // 위치
    double x = coeffs_x(0) + coeffs_x(1) * t + coeffs_x(2) * t * t + coeffs_x(3) * t * t * t;
    double y = coeffs_y(0) + coeffs_y(1) * t + coeffs_y(2) * t * t + coeffs_y(3) * t * t * t;

    // 1차 미분 (속도 벡터)
    double dx = coeffs_x(1) + 2 * coeffs_x(2) * t + 3 * coeffs_x(3) * t * t;
    double dy = coeffs_y(1) + 2 * coeffs_y(2) * t + 3 * coeffs_y(3) * t * t;

    // 2차 미분 (가속도 벡터)
    double ddx = 2 * coeffs_x(2) + 6 * coeffs_x(3) * t;
    double ddy = 2 * coeffs_y(2) + 6 * coeffs_y(3) * t;

    // 헤딩 ψ = atan2(dy, dx)
    double heading = std::atan2(dy, dx);

    // 곡률 κ = (dx * ddy - dy * ddx) / (dx² + dy²)^(3/2)
    double denom = std::pow(dx * dx + dy * dy, 1.5);
    double curvature = 0.0;
    if (denom > 1e-6) {
        curvature = (dx * ddy - dy * ddx) / denom;
    } else {
        curvature = 0.0;  // 속도가 너무 작으면 곡률 계산이 불안정
    }

    return SplineEval{Vector2d(x, y), heading, curvature};
}

bool isInsideTrack(double x, double y, const DMap& sampling_map) {

    if (sampling_map.empty() || sampling_map.begin()->second.empty()) {
        return false;
    }

    const auto& x_ref = sampling_map.at(__x_ref);            // 중심선 x
    const auto& y_ref = sampling_map.at(__y_ref);            // 중심선 y
    const auto& x_normvec = sampling_map.at(__x_normvec);    // 중심선에서의 노멀 x
    const auto& y_normvec = sampling_map.at(__y_normvec);    // 중심선에서의 노멀 y
    const auto& width_left = sampling_map.at(__width_left);  // 왼쪽 경계 폭
    const auto& width_right = sampling_map.at(__width_right);// 오른쪽 경계 폭

    double min_dist_sq = std::numeric_limits<double>::max();
    int closest_idx = -1;

     // 가장 가까운 기준선 인덱스 찾기
    for (size_t i = 0; i < x_ref.size(); ++i) {
        double dx = x - x_ref[i];
        double dy = y - y_ref[i];
        double dist_sq = dx * dx + dy * dy;
        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest_idx = i;
        }
    }

    if (closest_idx == -1) return false;

    // 기준선 정보
    double ref_x = x_ref[closest_idx];
    double ref_y = y_ref[closest_idx];
    double norm_x = x_normvec[closest_idx];
    double norm_y = y_normvec[closest_idx];
    double left = width_left[closest_idx];
    double right = width_right[closest_idx];

    // 횡방향 거리 계산 (노멀 벡터 방향)
    double lateral_offset = (x - ref_x) * norm_x + (y - ref_y) * norm_y;

    // 경계 안에 있는지 판별
    return (lateral_offset >= -left && lateral_offset <= right);
}

bool isValidSpline(const RowVector4d& coeffs_x,
                    const RowVector4d& coeffs_y,
                    double stepsize,
                    const DMap& sampling_map,
                    double kappa_max = 0.2) {


    const int num_samples = 20;

    for (int i = 0; i <= num_samples; ++i) {
        double t = stepsize * i / static_cast<double>(num_samples);
        SplineEval eval = evaluateSpline(coeffs_x, coeffs_y, t);

        // 트랙 안에 있는가?
        if (!isInsideTrack(eval.pos.x(), eval.pos.y(), sampling_map)) {
            return false;
        }

        // 곡률이 너무 큰가?
        if (std::abs(eval.curvature) > kappa_max) {
            return false;
        }        
    }

    return true;  // 모든 샘플이 유효하면 true
}

void genEdge(Graph& graph, 
    const NodeMap& nodesPerLayer, 
    const Offline_Params& params,
    const IVector& raceline_index_array,
    bool closed = true) {

        if (params.LAT_OFFSET <= 0.0) {
            cout << "Too small lateral offset" << endl;
        }

        const auto& xs = sampling_map[__x_raceline];
        const auto& ys = sampling_map[__y_raceline];

        Eigen::MatrixXd raceline_cl(xs.size(), 2);
        for (size_t i = 0; i < xs.size(); ++i) {
            raceline_cl(i, 0) = xs[i];
            raceline_cl(i, 1) = ys[i];
        }

        // raceline 전체를 스플라인으로 만들고 계수 가져오기
        auto raceline_result = calcSplines(raceline_cl);
        const MatrixXd& x_coeff_r = raceline_result.coeffs_x;
        const MatrixXd& y_coeff_r = raceline_result.coeffs_y;

        // closed -> 마지막과 0번 레이어 연결, open -> break
        for (size_t layer = 0; layer < nodesPerLayer.size(); ++layer) {
            
            const size_t start_layer = layer;
            const size_t end_layer = (layer + 1) % nodesPerLayer.size();
        
            if (!params.CLOSURE_DETECTION_DIST && end_layer == 0)
            break;

        // 1. 각 레이어의 레이스라인 기준 인덱스
        const int start_race_idx = raceline_index_array[start_layer];
        const int end_race_idx   = raceline_index_array[end_layer];

        // 2. 시작 레이어의 각 노드에 대해 반복
        const auto& start_layer_nodes = nodesPerLayer[start_layer];
        const auto& end_layer_nodes   = nodesPerLayer[end_layer];

        for (size_t start_idx = 0; start_idx < start_layer_nodes.size(); ++start_idx) {
            const Node& start_node = start_layer_nodes[start_idx];

        // 3. 레이스라인에서 같은 위치에 있는 end 레이어의 참조 노드 인덱스 계산
        int rel_race_offset = static_cast<int>(start_idx) - start_race_idx;
        int ref_end_idx = end_race_idx + rel_race_offset;

        // 4. 범위 클램프
        ref_end_idx = std::clamp(ref_end_idx, 0, static_cast<int>(end_layer_nodes.size() - 1));

        // 5. 거리 계산을 위한 노드 좌표 가져오기
        const Node& ref_end_node = end_layer_nodes[ref_end_idx];

        // 6. 거리 계산을 위한 좌표 행렬 생성
        MatrixXd spline_path(2, 2);
        spline_path << start_node.x, start_node.y,
                    ref_end_node.x, ref_end_node.y;


            VectorXd el_lengths = computeEuclideanDistances(spline_path);
            double dist = el_lengths(0);

            // 커브면 더 많이 연결
            double factor = (start_node.kappa > params.CURVE_THR) ? 2.0 : 1.0;
            int lat_steps = round(factor * dist * params.LAT_OFFSET / params.LAT_RESOLUTION);

            for (int destIdx = std::max(0, destIdx - lat_steps);
                destIdx <= std::min(static_cast<int>(nodesPerLayer[end_layer].size() - 1), destIdx + lat_steps);
                ++destIdx) {

                const Node& endNode = nodesPerLayer[end_layer][destIdx];
                MatrixXd x_coeffs, y_coeffs;

                if (start_node.raceline && endNode.raceline) {
                    x_coeffs = x_coeff_r.row(start_layer);
                    y_coeffs = y_coeff_r.row(start_layer);
                } else {
                    MatrixXd path(2, 2);
                    path << start_node.x, start_node.y,
                            endNode.x, endNode.y;

                    auto result = calcSplines(path, nullptr, start_node.psi, endNode.psi);

                    const MatrixXd& x_coeff_r = raceline_result.coeffs_x;
                    const MatrixXd& y_coeff_r = raceline_result.coeffs_y;
                }

                // 그래프에 엣지 추가
                ITuple src_key(start_layer, start_node.node_idx);
                graph.addEdge(src_key, endNode.node_idx); 
            }
        }  
    }
}

void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params) {
    plt::clf();

    // 트랙 경계선
    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // 레이싱 라인 및 샘플링된 포인트
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}, {"label", "Raceline"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}, {"label", "Sampled Raceline"}});
    plotHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);

    // 노드들
    plotHeading(nodesPerLayer);

    // --- Graph 엣지 (스플라인) 시각화 ---
    
    DVector spline_x_pts; 
    DVector spline_y_pts;

    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& current_node : layer_nodes) {
            ITuple src_key(current_node.layer_idx, current_node.node_idx);
            IVector child_nodes_idx;

            try {
                graph.getChildIdx(src_key, child_nodes_idx);
            } catch (const std::runtime_error& e) {
                continue;
            }
            
            for (int dest_node_idx : child_nodes_idx) {
                // 각 스플라인 그리기 전에 벡터를 비워줍니다. (이전 스플라인 점 데이터 초기화)
                spline_x_pts.clear(); 
                spline_y_pts.clear(); 
                
                size_t next_layer_idx = (current_node.layer_idx + 1) % nodesPerLayer.size();
                // next_node_idx 유효성 검사 추가 (인덱스 범위 체크)
                if (dest_node_idx < 0 || dest_node_idx >= nodesPerLayer[next_layer_idx].size()) {
                    // cout << "Warning: Invalid dest_node_idx " << dest_node_idx << " for layer " << next_layer_idx << endl; // 이 라인도 extended character 오류의 원인이 될 수 있습니다.
                    continue;
                }
                const Node& next_node = nodesPerLayer[next_layer_idx][dest_node_idx];

                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y,
                               next_node.x, next_node.y;
                
                double psi_s = current_node.psi;
                double psi_e = next_node.psi;

                VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

                SplineResult res;
                try {
                    // calcSplines의 마지막 인자인 use_dist_scaling은 true로 가정
                    res = calcSplines(spline_path, &el_lengths, psi_s, psi_e, true);
                } catch (const std::exception& e) {
                    continue;
                }

                // 계산된 스플라인 계수를 사용하여 곡선 그리기
                const int num_spline_segments = 10; 
                for (int k = 0; k <= num_spline_segments; ++k) {
                    double t_eval = static_cast<double>(k) / num_spline_segments;
                    SplinePoint sp = evaluateSpline(res.coeffs_x.row(0), res.coeffs_y.row(0), t_eval, res.ds(0), true);
                    spline_x_pts.push_back(sp.x);
                    spline_y_pts.push_back(sp.y);
                }
                plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}}); // {"label", "Valid Splines"}
            }
        }
    }

    plt::title("Track and Planned Graph");
    plt::grid(true);
    plt::axis("equal");
    plt::legend();
    plt::show();
}



int main() {
    IVector idx_sampling;
    Offline_Params params;

    string map_file_in = "inputs/gtpl_levine.csv";
    string map_file_out = "inputs/gtpl_levine_out.csv";

    // global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    readDMapFromCSV(map_file_in, gtpl_map);

    addDVectorToMap(gtpl_map, "bound_r");
    addDVectorToMap(gtpl_map, "bound_l");
    addDVectorToMap(gtpl_map, "raceline");
    addDVectorToMap(gtpl_map, "delta_s");

    writeDMapToCSV(map_file_out, gtpl_map);
    
    // layer 간격을 위한 raceline points sampling 
    samplePointsFromRaceline(gtpl_map[__kappa],
                             gtpl_map[__delta_s],
                             params.LON_CURVE_STEP,
                             params.LON_STRAIGHT_STEP,
                             params.CURVE_THR,
                             idx_sampling);

    // cout << "idx size:" << idx_sampling.size() << endl;
    
    for (const auto& [key, vec] : gtpl_map) {
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < vec.size()) {
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }
    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);

    calcHeading(sampling_map[__x_raceline],
                sampling_map[__y_raceline],
                sampling_map[__psi]);
    // 여기서 계산되는 sampling된 bound_l, r은 node 생성 시에만 쓰인다. 
    calcHeading(sampling_map[__x_bound_l],
                sampling_map[__y_bound_l],
                sampling_map[__psi_bound_l]);
    calcHeading(sampling_map[__x_bound_r],
                sampling_map[__y_bound_r],
                sampling_map[__psi_bound_r]);  

    // 노드 그리드 생성
    NodeMap nodesPerLayer;
    IVector raceline_index_array;
    Vector2d node_pos;

    genNode(nodesPerLayer,
            raceline_index_array,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    //Graph edgeList;

    MatrixXd path_xy(sampling_map[__x_raceline].size(), 2);
    for (size_t i = 0; i < sampling_map[__x_raceline].size(); ++i) {
        path_xy(i, 0) = sampling_map[__x_raceline][i];
        path_xy(i, 1) = sampling_map[__y_raceline][i];
    }

    double psi_s = sampling_map[__psi][0];
    double psi_e = sampling_map[__psi].back();

    SplineResult result = calcSplines(path_xy, nullptr, psi_s, psi_e);

    // 최종 그래프 생성
    Graph myGraph;
    genEdge(myGraph,
            nodesPerLayer,
            params,
            raceline_index_array);
    
    // 시각화
    // visual(nodesPerLayer, myGraph, params);

    return 0;
}