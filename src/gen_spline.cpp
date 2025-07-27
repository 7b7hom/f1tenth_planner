#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// --- 실행하기 위해서 main에서 가져온 부분 ---

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
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    size_t i = 0;
    for (const auto& [key, _] : map) {
        file << key;
        if (++i != num_cols) file << delimiter;
    }
    file << '\n';
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

// DVector를 Map 구조로 추가(연산)
void addDVectorToMap(DMap &map, string attr, const IVector *idx_array = nullptr) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();
    } else {
        len = idx_array->size();
    }
    DVector x_out(len), y_out(len);
    string x_label = "x_" + attr;
    string y_label = "y_" + attr;
    
    if (!attr.compare("bound_r")) {
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__width_right][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__width_right][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    } else if (!attr.compare("bound_l")) {
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] - map[__x_normvec][i] * map[__width_left][i];
            y_out[i] = map[__y_ref][i] - map[__y_normvec][i] * map[__width_left][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    } else if (!attr.compare("raceline")) {
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__alpha][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__alpha][i];
        }
        map[x_label] = x_out;
        map[y_label] = y_out;
    } else if (!attr.compare("delta_s")) {
        for (size_t i = 0; i < len - 1; ++i) {
            x_out[i] = map[__s_racetraj][i+1] - map[__s_racetraj][i];
        }
        map[attr] = x_out;
    }
}

// 레이싱 라인에서 경로 계획을 위한 layer가 될 지점들을 샘플링
void samplePointsFromRaceline(const DVector& kappa, const DVector& dist,
                              double d_curve, double d_straight, double curve_th, IVector& idx_array) {
    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist;
        }
        if ((cur_dist + dist[i]) > next_dist) {
            idx_array.push_back(static_cast<int>(i));
            if (fabs(kappa[i]) < curve_th) {
                next_dist += d_straight;
            } else {
                next_dist += d_curve;
            }
            next_dist_min = cur_dist + d_curve;
        }
        cur_dist += dist[i];
    }
}

// 주어진 각도를 -PI에서 PI 사이의 값으로 정규화
double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// 주어진 X, Y 좌표 벡터를 기반으로 각 지점에서의 헤딩(Heading, 진행 방향 각도)을 계산
void calcHeading(DVector &x_raceline, DVector &y_raceline, DVector &psi) {
    size_t N = x_raceline.size();
    psi.resize(N);
    double dx, dy;
    for (size_t i = 0; i < N; ++i) {
        if (i != N -1) {
            dx = x_raceline[i+1] - x_raceline[i];
            dy = y_raceline[i+1] - y_raceline[i];
        } else { // 닫힌 회로 가정
            dx = x_raceline[0] - x_raceline[N - 1];
            dy = y_raceline[0] - y_raceline[N - 1];
        }
        psi[i] = atan2(dy, dx) - M_PI_2;
        normalizeAngle(psi[i]);
    }
}

// 샘플링된 레이어마다 경로 계획을 위한 Node(차량이 횡방향으로 이동 가능한 위치들) 생성
void genNode(NodeMap& nodesPerLayer, const double veh_width, float lat_resolution) {
    const size_t N = sampling_map[__alpha].size();
    nodesPerLayer.resize(N); 

    for (size_t i = 0; i < N; ++i){ // 각 레이어(층)에 대해 반복

        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        
        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;
        nodesPerLayer[i].resize(num_nodes); // 현재 레이어의 노드 벡터 크기 조정 (NodeMap은 vector<vector<Node>> 이므로 inner vector의 resize)
        // cout << i << "번째 레이어의 노드 개수: " << num_nodes << endl;
        for (int idx = 0; idx < num_nodes; ++idx) { 
            double alpha = start_alpha + idx * lat_resolution; // 현재 노드의 횡방향 오프셋 계산
            Vector2d node_pos = ref_xy + alpha * norm_vec; // 노드의 (x, y) 좌표 계산 (여기서 선언)

            Node current_node_instance; // 각 노드 인스턴스를 이 루프 안에서 새로 생성하여 초기화 문제를 방지
            
            // 필수 멤버 초기화 및 할당
            current_node_instance.layer_idx = i; // 현재 레이어 인덱스
            current_node_instance.node_idx = idx; // 노드 인덱스
            current_node_instance.x = node_pos.x();
            current_node_instance.y = node_pos.y();
            current_node_instance.kappa = 0.0;
            current_node_instance.raceline = (idx == raceline_index); // raceline_index는 현재 레이어의 레이싱 라인 인덱스


            // --- 노드의 헤딩(psi) 계산 (보간) ---
            double psi_interp;
            if (idx < raceline_index) { 
                if (abs(sampling_map[__psi_bound_l][i] - sampling_map[__psi][i]) >= M_PI) {
                    double bl = sampling_map[__psi_bound_l][i] + 2 * M_PI * (sampling_map[__psi_bound_l][i] < 0);
                    double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
                    psi_interp = bl + (p - bl) * idx / raceline_index; 
                } else {
                    psi_interp = sampling_map[__psi_bound_l][i] + (sampling_map[__psi][i] - sampling_map[__psi_bound_l][i]) * (idx+1) / raceline_index;
                }
                current_node_instance.psi = normalizeAngle(psi_interp);
            }
            else if (idx == raceline_index) { 
                psi_interp = sampling_map[__psi][i];
                current_node_instance.psi = psi_interp;
            }
            else { 
                int remain = num_nodes - raceline_index - 1;
                double t = static_cast<double>(idx - raceline_index) / std::max(remain, 1); 
                psi_interp = sampling_map[__psi][i] + t * (sampling_map[__psi_bound_r][i] - sampling_map[__psi][i]);
                current_node_instance.psi = normalizeAngle(psi_interp);
            }
            //current_node_instance.psi = sampling_map[__psi][i];
            
            nodesPerLayer[i][idx] = current_node_instance; // <-- 생성된 노드 인스턴스를 NodeMap에 할당
        }
    }
}


// 주어진 X, Y 좌표에서 해당 psi (헤딩) 방향을 화살표로 시각화
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
        plt::plot(x_line, y_line, {{"color", "red"}});

        #if 1 // 화살촉 그리기
        theta = atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);
        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "red"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "red"}});
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


// --- spline 관련 함수들 ---
// 여기부터 코딩

VectorXd computeEuclideanDistances(const MatrixXd& path) {
    int N = path.rows() - 1;
    VectorXd dists(N);
    for(int i = 0; i < N; ++i){
        dists(i) = (path.row(i + 1) - path.row(i)).norm();
    }
    return dists;
}

SplineResult calcSplines(const MatrixXd& path, // spline 생성 시 기준이 되는 경로 점들의 X, Y 좌표 담고 있는 참조 변수(&: 오버헤드 줄여줌.)
                         const VectorXd* el_lengths_ptr = nullptr, // 유클리드 거리 담고 있는 포인터 변수(nullptr이면 유클리드 거리 직접 계산, 아니면 포인터가 가리키는 VectorXd 객체 거리 사용)
                         double psi_s = std::numeric_limits<double>::quiet_NaN(), // spline 시작점 heading(NaN이면 헤딩 지정x -> natural spline 조건 따름(2차 미분값 0))
                         double psi_e = std::numeric_limits<double>::quiet_NaN(), // spline 끝점 heading
                         bool use_dist_scaling = true){ // spline의 1차 및 2차 미분 연속성 조건에 거리 스케일링 적용할지 여부(기본값 true -> el_lengths 고려하여 미분값들 스케일링)
    
    bool closed = (path.row(0) - path.bottomRows(1)).norm() < 1e-6;

    int no_splines = path.rows() - 1;

    VectorXd el_lengths;
    if(use_dist_scaling && el_lengths_ptr == nullptr){
        el_lengths = computeEuclideanDistances(path);
    }else if(el_lengths_ptr){

    }

    if(use_dist_scaling && closed){
        VectorXd tmp(el_lengths.size() + 1);
        tmp << el_lengths, el_lengths(0);
        el_lengths = tmp;
    }

    VectorXd scaling = VectorXd::Ones(no_splines - 1);
    if(use_dist_scaling && no_splines > 1){
        scaling = el_lengths.head(no_splines - 1).array() / el_lengths.segment(1, no_splines - 1).array();
    }

    const int dim = no_splines * 4; 
    MatrixXd M = MatrixXd::Zero(dim, dim);
    VectorXd b_x = VectorXd::Zero(dim);
    VectorXd b_y = VectorXd::Zero(dim);

    Matrix<double, 4, 8> template_M;
    template_M << 1, 0, 0, 0, 0, 0, 0, 0,
                  1, 1, 1, 1, 0, 0, 0, 0, 
                  0, 1, 2, 3, 0, -1, 0, 0,
                  0, 0, 2, 6, 0, 0, -2, 0;
    
    for(int i = 0; i < no_splines; ++i){
        int j= i * 4;

        if(i < no_splines - 1){
            M.block(j, j, 4, 8) = template_M;
            M(j + 2, j + 5) *= scaling(i);
            M(j + 3, j + 6) *= pow(scaling(i), 2);
        }else{
            M.block(j, j, 2, 4) << 1, 0, 0, 0,
                                   1, 1, 1, 1;
        }

        b_x.segment(j, 2) << path(i, 0), path(i + 1, 0);
        b_y.segment(j, 2) << path(i, 1), path(i + 1, 1);
    }

    if(!closed){
        double el_length_s = el_lengths_ptr ? (*el_lengths_ptr)(0) : 1.0;
        double el_length_e = el_lengths_ptr ? el_lengths_ptr->tail(1)(0) : 1.0;

        M(dim - 2, 1) = 1.0;
        b_x(dim - 2) = cos(psi_s + M_PI_2) * el_length_s;
        b_y(dim - 2) = sin(psi_s + M_PI_2) * el_length_s;

        M.block(dim - 1, dim - 4, 1, 4) << 0, 1, 2, 3;
        b_x(dim - 1) = cos(psi_e + M_PI_2) * el_length_e;
        b_y(dim - 1) = sin(psi_e + M_PI_2) * el_length_e;
    }else{
        M(dim - 2, 1) = scaling.tail(1)(0);
        M.block(dim - 2, dim - 3, 1, 3) << -1, -2, -3;
        M(dim - 1, 2) = 2 * pow(scaling.tail(1)(0), 2);
        M.block(dim - 1, dim - 2, 1, 2) << -2, -6;
    }

    VectorXd x_les = M.fullPivLu().solve(b_x);
    VectorXd y_les = M.fullPivLu().solve(b_y);

    MatrixXd coeffs_x(no_splines, 4), coeffs_y(no_splines, 4);
    for (int i = 0; i < no_splines; ++i) {
        coeffs_x.row(i) = x_les.segment(i * 4, 4).transpose();
        coeffs_y.row(i) = y_les.segment(i * 4, 4).transpose();
    }

    MatrixXd normvec(no_splines, 2);
    for (int i = 0; i < no_splines; ++i) {
        double dx = coeffs_y(i, 1);
        double dy = -coeffs_x(i, 1);
        double norm = std::sqrt(dx * dx + dy * dy);
        normvec(i, 0) = dx / norm;
        normvec(i, 1) = dy / norm;
    }

    VectorXd ds(no_splines);
    for (int i = 0; i < no_splines; ++i)
        ds(i) = (path.row(i + 1) - path.row(i)).norm();

    return SplineResult{
        .coeffs_x = coeffs_x,
        .coeffs_y = coeffs_y,
        .M = M,
        .normvec_normalized = normvec,
        .ds = ds
    };

}

// spline 평가하여 위치, 헤딩, 곡률 반환
SplinePoint evaluateSpline(const RowVector4d& coeff_x, const RowVector4d& coeff_y,
                           double t, double ds_current, bool normalized_t){
    SplinePoint sp;

    // spline 계수들은 t가 0~1로 가정하고 계산하므로 정규화 필요
    double local_t = t;
    if(!normalized_t){
        local_t = t / ds_current;
    }

    // 위치 계산
    sp.x = coeff_x(0) + coeff_x(1) * local_t + coeff_x(2) * pow(local_t, 2) + coeff_x(3) * pow(local_t, 3);
    sp.y = coeff_y(0) + coeff_y(1) * local_t + coeff_y(2) * pow(local_t, 2) + coeff_y(3) * pow(local_t, 3);
    
    // 1차 미분
    double dx_dt = coeff_x(1) + 2 * coeff_x(2) * local_t + 3 * coeff_x(3) * pow(local_t, 2);
    double dy_dt = coeff_y(1) + 2 * coeff_y(2) * local_t + 3 * coeff_y(3) * pow(local_t, 2);
    
    // heading 계산
    sp.psi = normalizeAngle(atan2(dy_dt, dx_dt));

    // 2차 미분
    double d2x_dt2 = 2 * coeff_x(2) + 6 * coeff_x(3) * local_t;
    double d2y_dt2 = 2 * coeff_y(2) + 6 * coeff_y(3) * local_t;

    // curvature 계산 (2D 매개변수 곡선의 곡률 공식 이용)
    // 추가적인 명시적 스케일링 팩터(1/ds_current, 1/ds_current^2)가 필요 x(앞에서 local_t = t/ds_current로 이미 내포)
    double numerator = dx_dt * d2y_dt2 - dy_dt * d2x_dt2;
    double denominator = pow(pow(dx_dt, 2) + pow(dy_dt, 2), 1.5);
    if(denominator < 1e-9){  // 0에 가까울 때 0으로 나누는 오류 방지
        sp.kappa = 0.0;
    }else{
        sp.kappa = numerator / denominator;
    }

    return sp;
}

// 주어진 (x, y) 점이 트랙 경계 내부에 있는지(왼, 오 경계 사이에 위치하는지) 확인하는 함수 (전역 sampling_map 사용)
bool isPointInsideTrackBounds(double x, double y){
    if(sampling_map.empty() || sampling_map.begin()->second.empty()){
        return false;
    }

    // (x, y) 점에서 가장 가까운 sampling_map의 기준선 인덱스 찾기
    double min_dist_sq = std::numeric_limits<double>::max(); // 가능한 가장 큰 수로 초기갑 설정
    int closest_ref_idx = -1; // 유효하지 않은 값으로 초기값 설정

    for(size_t i = 0; i < sampling_map[__x_ref].size(); ++i){
        double dx = x - sampling_map[__x_ref][i];
        double dy = y - sampling_map[__y_ref][i];
        double dist_sq = dx * dx + dy * dy; // sprt() 연산 피하기(성능 위해)
        if(dist_sq < min_dist_sq){
            min_dist_sq = dist_sq;
            closest_ref_idx = i;
        }
    }

    if(closest_ref_idx == -1){
        return false;
    }

    // 가장 가까운 기준선 인덱스의 정보 가져오기
    double ref_x = sampling_map[__x_ref][closest_ref_idx]; // __x_ref, __y_ref: 기준선 좌표
    double ref_y = sampling_map[__y_ref][closest_ref_idx];
    double norm_x = sampling_map[__x_normvec][closest_ref_idx]; // __x_normvec, __y_normvec: 기준선 법선 벡터
    double norm_y = sampling_map[__y_normvec][closest_ref_idx];
    double width_left = sampling_map[__width_left][closest_ref_idx]; // __width_left, __width_right: 기준선으로부터의 좌우 폭
    double width_right = sampling_map[__width_right][closest_ref_idx];

    // spline 점의 기준선 법선 방향 횡방향 오프셋 계산
    // -> (x, y) 점이 기준선으로부터 법선 벡터 방향으로 얼마나 떨어져 있는가
    double lateral_offset = (x - ref_x) * norm_x + (y - ref_y) * norm_y; // 법선 벡터 norm_x, norm_y: 기준선에 수직인 방향을 가리키는 단위 벡터
    /*cout << "DEBUG OOB: Point(" << x << "," << y << ") RefIdx=" << closest_ref_idx
     << " Offset=" << lateral_offset << " Bounds=[" << -width_left << "," << width_right << "]"
     << " NormVec=(" << norm_x << "," << norm_y << ")" << endl;*/

    double epsilon = 1e-6;
    // lateral_offset이 트랙의 유효한 횡방향 범위 내에 있는지 직접적으로 확인하는 최종 단계
    if(lateral_offset >= -width_left - epsilon && lateral_offset <= width_right + epsilon){
        return true;
    }else{
        return false;
    }
}

// 하나의 스플라인 구간이 유효한 경로로 사용될 수 있는지 검사
bool checkSplineValidity(const RowVector4d coeff_x, const RowVector4d& coeff_y, double ds_current,
                         const Offline_Params& params){ // const MatrixXd& track_bounds
    // spline 경로 샘플링
    const int num_samples = 10;

    double max_allowed_kappa = 30.0 / params.VEH_TURN;

    for(int k = 0; k <= num_samples; ++k){
        double t_eval = static_cast<double>(k) / num_samples; // 파라미터 t 값을 균등하게 분할
        SplinePoint sp = evaluateSpline(coeff_x, coeff_y, t_eval, ds_current, true);

        // 트랙 경계 벗어나는지 확인
        if(!isPointInsideTrackBounds(sp.x, sp.y)){
            cout << "REJECTED (OUT OF BOUNDS): Spline from " << sp.x << "," << sp.y << " is out of bounds." << endl;
            return false;
        }

        // 곡률 제약 조건 확인
        // std::abs(sp.kappa): 현재 spline 점에서의 kappa 절댓값
        // 4.0 / params.VEH_TURN: 차량이 허용하는 최대 곡률 (최소 회전 반경의 역수)
        if(std::abs(sp.kappa) > max_allowed_kappa){
            cout << "REJECTED (EXCESSIVE CURVATURE): Point (" << sp.x << "," << sp.y << "), kappa=" << sp.kappa << ", max_allowed=" << max_allowed_kappa << endl;
            return false;
        }
    }

    if (true) { // 모든 스플라인에 대해 kappa 확인용 로그
    double max_kappa = 0.0;
    double min_kappa = 1e9;
    for (int k = 0; k <= num_samples; ++k) {
        double t_eval = static_cast<double>(k) / num_samples;
        SplinePoint sp = evaluateSpline(coeff_x, coeff_y, t_eval, ds_current, true);
        max_kappa = std::max(max_kappa, std::abs(sp.kappa));
        min_kappa = std::min(min_kappa, std::abs(sp.kappa));
    }
    std::cout << "[PASS] SPLINE OK: max_kappa=" << max_kappa
              << ", min_kappa=" << min_kappa
              << ", limit=" << max_allowed_kappa
              << std::endl;
    }

    return true;
}

void generateGraphEdges(Graph& graph, const NodeMap& nodesPerLayer, const Offline_Params params){
    const size_t num_layers = nodesPerLayer.size();

    // layer 순회
    for(size_t current_layer_idx = 0; current_layer_idx < num_layers; ++current_layer_idx){
        size_t next_layer_idx = (current_layer_idx + 1) % num_layers;

        const auto& current_nodes_in_layer = nodesPerLayer[current_layer_idx];
        const auto& next_nodes_in_layer = nodesPerLayer[next_layer_idx];

        // 현재 layer의 node 순회
        for(const auto& current_node : current_nodes_in_layer){
            // lat_steps 로직
            int refDestIdx = current_node.node_idx; // 기준 목적지 인덱스 선정(current_node와 같은 횡방향 인덱스를 가진 다음 layer의 node)
            refDestIdx = clamp(refDestIdx, 0, static_cast<int>(next_nodes_in_layer.size() - 1));
            
            const Node& refEndNode = next_nodes_in_layer[refDestIdx];

            // 현재 layer와 다음 layer 기준 노드 사이의 거리 계산
            double dist_between_layers = (Vector2d(refEndNode.x, refEndNode.y) - Vector2d(current_node.x, current_node.y)).norm();

            double ratio = 0.0;
            if(params.CURVE_THR > 1e-9){
                ratio = min(abs(current_node.kappa) / params.CURVE_THR, 2.0); // 현재 노드의 곡률을 CURVE_THR로 정규화
            }
            // factor -> 직선보다 곡선에서 더 넓은 노드 탐색을 가능하게 하는 계수
            double factor = 1.0 / (1.0 + 0.5 * ratio); // 곡률 높을수록 lat_steps 줄이기

            // 곡률이 높고 거리가 멀수록 더 많은 노드를 살펴봄
            int lat_steps = static_cast<int>(round(factor * dist_between_layers * params.LAT_OFFSET / params.LAT_RESOLUTION));

            lat_steps = min(lat_steps, (int)params.MAX_LAT_STEPS);

            for(int dest_node_idx = max(0, refDestIdx - lat_steps); dest_node_idx <= min(static_cast<int>(next_nodes_in_layer.size() - 1), refDestIdx + lat_steps); ++dest_node_idx){
                const Node& next_node = next_nodes_in_layer[dest_node_idx];

                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y, next_node.x, next_node.y;

                double psi_s = current_node.psi;
                double psi_e = next_node.psi;

                VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

                SplineResult res;
                try{
                    res = calcSplines(spline_path, &el_lengths, psi_s, psi_e, true);
                }catch(const exception& e){
                    continue;
                }

                if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params)){
                    ITuple src_key(current_node.layer_idx, current_node.node_idx);
                    graph.addEdge(src_key, next_node.node_idx);
                    cout << "SPLINE PASSED!!!! from (" << current_node.layer_idx << "," << current_node.node_idx
                         << ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }else{
                    cout << "SPLINE REJECTED from (" << current_node.layer_idx << "," << current_node.node_idx
                         << ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }
            }

        }
    }
}


// -- prune_graph() --

void prune_graph(Graph& graph, int num_layers, bool closed = true){
    int j = 0;
    int rmv_cnt_tot = 0;

    vector<ITuple> nodes;
    for(int layer = 0; layer < num_layers; ++layer){
        for(const auto& [key, vec] : graph.getAdjLists()){
            if(get<0>(key) == layer){
                nodes.push_back(key);
            }
        }
    }

    while(true){
        int rmv_cnt = 0;

        for(const auto& node : nodes){
            int layer = get<0>(node);
            int node_idx = get<1>(node);

            if(!closed && (layer == 0 || layer == num_layers - 1)){
                continue;
            }

            IVector children;
            vector<ITuple> parents;

            try{
                graph.getChildIdx(node, children);
            }catch(...){
                children.clear();
            }
            
            graph.getParentNode(layer, node_idx, parents);

            if(children.empty() || parents.empty()){
                if(children.empty()){
                    for(const auto& parent : parents){
                        graph.removeEdge(parent, node_idx);
                        ++rmv_cnt;
                    }
                }else{
                    for(int child : children){
                        graph.removeEdge(node, child);
                        ++rmv_cnt;
                    }
                }
            }
        }

        if(rmv_cnt == 0){
            break;
        }else{
            rmv_cnt_tot += rmv_cnt;
        }
        ++j;
    }

}

// 트랙의 경계, 레이싱 라인, 샘플링된 포인트, 생성된 노드들, 그리고 그래프 엣지(스플라인)를 시각화
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

                plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}});
                /*if(current_node.layer_idx==0){
                    if(dest_node_idx==0){
                        plt::plot(spline_x_pts, spline_y_pts, {{"color", "yellow"}, {"linewidth", "1"}}); // {"label", "Valid Splines"}
                    }else{
                        plt::plot(spline_x_pts, spline_y_pts, {{"color", "black"}, {"linewidth", "1"}}); // {"label", "Valid Splines"}
                    }
                }else{
                    plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}}); // {"label", "Valid Splines"}
                }*/
            }
        }
    }

    plt::title("Track and Planned Graph");
    plt::grid(true);
    plt::axis("equal");
    plt::legend();
    plt::show();
}

// 전체 경로 계획 파이프라인을 실행하는 함수
void runPlanningPipeline(const Offline_Params& params, const std::string& map_file_in, const std::string& map_file_out) {
    // 1. 트랙 데이터 로드 및 전처리
    readDMapFromCSV(map_file_in, gtpl_map); // gen_spline.cpp의 전역 gtpl_map에 로드
    addDVectorToMap(gtpl_map, "bound_r");
    addDVectorToMap(gtpl_map, "bound_l");
    addDVectorToMap(gtpl_map, "raceline");
    addDVectorToMap(gtpl_map, "delta_s");
    writeDMapToCSV(map_file_out, gtpl_map);

    // 2. 레이어 샘플링
    IVector idx_sampling;
    samplePointsFromRaceline(gtpl_map[__kappa], gtpl_map[__delta_s],
                             params.LON_CURVE_STEP, params.LON_STRAIGHT_STEP,
                             params.CURVE_THR, idx_sampling);
    
    // 샘플링된 인덱스를 사용하여 gtpl_map에서 데이터를 복사하여 sampling_map 채우기
    for (const auto& [key, vec] : gtpl_map) {
        sampling_map[key].reserve(idx_sampling.size()); 
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < (int)vec.size()) { 
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }
    addDVectorToMap(sampling_map, "delta_s", &idx_sampling); 
    
    calcHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);
    calcHeading(sampling_map[__x_bound_l], sampling_map[__y_bound_l], sampling_map[__psi_bound_l]);
    calcHeading(sampling_map[__x_bound_r], sampling_map[__y_bound_r], sampling_map[__psi_bound_r]);

    // 3. 노드 그리드 생성
    NodeMap nodesPerLayer;
    genNode(nodesPerLayer, params.VEH_WIDTH, params.LAT_RESOLUTION);

    // 4. 스플라인 생성 및 유효성 검사, 최종 그래프 구축
    Graph directedGraph; 
    generateGraphEdges(directedGraph, nodesPerLayer, params);

    prune_graph(directedGraph, nodesPerLayer.size(), true); 

    // 5. 최종 그래프 연결 확인 (Print Graph)
    cout << "\n--- 최종 생성된 그래프 (유효한 스플라인 엣지 포함) ---" << endl;
    directedGraph.printGraph();

    // 6. 결과 시각화 (Plotting)
    visual(nodesPerLayer, directedGraph, params); 
}

