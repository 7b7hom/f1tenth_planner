#include "graph_planner.hpp"
#include "config_modena.h"

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

// DVector를 Map 구조로 추가
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
        psi[i] = normalizeAngle(psi[i]);
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
                double t = static_cast<double>(idx - raceline_index) / max(remain, 1); 
                psi_interp = sampling_map[__psi][i] + t * (sampling_map[__psi_bound_r][i] - sampling_map[__psi][i]);
                current_node_instance.psi = normalizeAngle(psi_interp);
            }
            //current_node_instance.psi = sampling_map[__psi][i];
            
            nodesPerLayer[i][idx] = current_node_instance; // <-- 생성된 노드 인스턴스를 NodeMap에 할당
        }
    }
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
                         double psi_s = numeric_limits<double>::quiet_NaN(), // spline 시작점 heading(NaN이면 헤딩 지정x -> natural spline 조건 따름(2차 미분값 0))
                         double psi_e = numeric_limits<double>::quiet_NaN(), // spline 끝점 heading
                         bool use_dist_scaling = true){ // spline의 1차 및 2차 미분 연속성 조건에 거리 스케일링 적용할지 여부(기본값 true -> el_lengths 고려하여 미분값들 스케일링)
    
    bool closed = (path.row(0) - path.bottomRows(1)).norm() < 1e-6;

    int no_splines = path.rows() - 1;

    VectorXd el_lengths;
    if(use_dist_scaling && el_lengths_ptr == nullptr){
        el_lengths = computeEuclideanDistances(path);
    }else if(el_lengths_ptr){
        el_lengths = *el_lengths_ptr;
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
        double norm = sqrt(dx * dx + dy * dy);
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
    double min_dist_sq = numeric_limits<double>::max(); // 가능한 가장 큰 수로 초기갑 설정
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

    double max_allowed_kappa = 30.0 / params.veh_turn;

    for(int k = 0; k <= num_samples; ++k){
        double t_eval = static_cast<double>(k) / num_samples; // 파라미터 t 값을 균등하게 분할
        SplinePoint sp = evaluateSpline(coeff_x, coeff_y, t_eval, ds_current, true);

        // 트랙 경계 벗어나는지 확인
        if(!isPointInsideTrackBounds(sp.x, sp.y)){
            cout << "REJECTED (OUT OF BOUNDS): Spline from " << sp.x << "," << sp.y << " is out of bounds." << endl;
            return false;
        }

        // 곡률 제약 조건 확인
        // abs(sp.kappa): 현재 spline 점에서의 kappa 절댓값
        // 4.0 / params.VEH_TURN: 차량이 허용하는 최대 곡률 (최소 회전 반경의 역수)
        if(abs(sp.kappa) > max_allowed_kappa){
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
        max_kappa = max(max_kappa, abs(sp.kappa));
        min_kappa = min(min_kappa, abs(sp.kappa));
    }
    /*cout << "[PASS] SPLINE OK: max_kappa=" << max_kappa
              << ", min_kappa=" << min_kappa
              << ", limit=" << max_allowed_kappa
              << endl;*/
    }

    return true;
}

void genEdges(Graph& graph, const NodeMap& nodesPerLayer, const Offline_Params params){
    const size_t num_layers = nodesPerLayer.size();

    auto& adj = graph.getAdjLists_mutable();
    for (size_t l = 0; l < num_layers; ++l) {
        for (size_t j = 0; j < nodesPerLayer[l].size(); ++j) {
            adj[ ITuple((int)l, (int)j) ]; // touch: 빈 map 생성
        }
    }

    // raceline 추종 엣지 살리기
    for(size_t current_layer_idx = 0; current_layer_idx < num_layers; ++current_layer_idx){
        size_t next_layer_idx = (current_layer_idx + 1) % num_layers;

        // 현재 layer의 레이스 라인 노드 찾기
        const Node* current_raceline_node = nullptr;
        for(const auto& node : nodesPerLayer[current_layer_idx]){
            if(node.raceline){
                current_raceline_node = &node;
                break;
            }
        }

        // 다음 layer의 레이스 라인 노드 찾기
        const Node* next_raceline_node = nullptr;
        for(const auto& node : nodesPerLayer[next_layer_idx]){
            if(node.raceline){
                next_raceline_node = &node;
                break;
            }
        }

        // 두 레이어 모두 레이스 라인 노드가 존재하면 edge 생성 시도
        if(current_raceline_node && next_raceline_node){
            MatrixXd spline_path(2, 2);
            spline_path << current_raceline_node->x, current_raceline_node->y, next_raceline_node->x, next_raceline_node->y;
            VectorXd el_lengths(1);
            el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

            try{
                SplineResult res = calcSplines(spline_path, &el_lengths, current_raceline_node->psi, next_raceline_node->psi, true);
                if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params)){
                    ITuple src_key(current_raceline_node->layer_idx, current_raceline_node->node_idx);
                    graph.addEdge(src_key, next_raceline_node->node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
                }                
            }catch(const exception& e){}
        }
    }

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
            if(params.curve_thr > 1e-9){
                ratio = min(abs(current_node.kappa) / params.curve_thr, 2.0); // 현재 노드의 곡률을 CURVE_THR로 정규화
            }
            // factor -> 직선보다 곡선에서 더 넓은 노드 탐색을 가능하게 하는 계수
            double factor = 1.0 / (1.0 + 0.5 * ratio); // 곡률 높을수록 lat_steps 줄이기

            // 곡률이 높고 거리가 멀수록 더 많은 노드를 살펴봄
            int lat_steps = static_cast<int>(round(factor * dist_between_layers * params.lat_offset / params.lat_resolution));

            lat_steps = min(lat_steps, (int)params.max_lat_steps);

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
                    graph.addEdge(src_key, next_node.node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0)); 
                    //cout << "SPLINE PASSED!!!! from (" << current_node.layer_idx << "," << current_node.node_idx
                         //<< ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }else{
                    //cout << "SPLINE REJECTED from (" << current_node.layer_idx << "," << current_node.node_idx
                         //<< ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }
            }

        }
    }
}

void prune_graph(Graph& graph, int num_layers, bool closed = true) {
    int total_removed_edges = 0;
    int iteration_count = 0;

    while (true) {
        int removed_edges_in_iter = 0;
        set<pair<ITuple, int>> edges_to_remove;
        vector<ITuple> all_current_graph_keys;

        for (const auto& [key, _] : graph.getAdjLists()) {
            all_current_graph_keys.push_back(key);
        }
        for (const auto& node_key : all_current_graph_keys) {
            int layer = get<0>(node_key);
            int node_idx = get<1>(node_key);

            if (!closed && (layer == 0 || layer == num_layers - 1)) {
                continue;
            }

            IVector children_idx;
            vector<ITuple> parents;

            try {
                graph.getChildIdx(node_key, children_idx);
            } catch (...) {
                children_idx.clear();
            }

            graph.getParentNode(layer, node_idx, parents, num_layers);

            if (parents.empty()) {
                for (int child_idx : children_idx) {
                    edges_to_remove.insert({node_key, child_idx});
                }
            }
            else if (children_idx.empty()) {
                for (const auto& parent : parents) {
                    edges_to_remove.insert({parent, node_idx});
                }
            }
        }

        for (const auto& edge : edges_to_remove) {
            graph.removeEdge(edge.first, edge.second);
            removed_edges_in_iter++;
        }

        cout << "Iteration " << iteration_count << ": Removed " << removed_edges_in_iter << " edges." << endl;

        if (removed_edges_in_iter == 0) {
            break;
        }

        total_removed_edges += removed_edges_in_iter;
        iteration_count++;
    }

    cout << "Total edges removed during pruning: " << total_removed_edges << endl;
}


void gen_offline_cost(Graph& graph, const Offline_Params& params, const NodeMap& nodesPerLayer){
    int total_edges_count = 0;
    for(const auto& [srcKey, destMap] : graph.getAdjLists()){
        total_edges_count += destMap.size();
    }
    
    int processed_edges_count = 0; 

    for(auto& [srcKey, destMap] : graph.getAdjLists_mutable()){ // mutable 접근자 통해 EdgeInfo 수정
        for(auto& [destIdx, edgeInfo] : destMap){ // EdgeInfo를 참조로 받아 직접 수정
            double current_offline_cost = 0.0; 

            // spline 곡률 계산을 위한 샘플링
            vector<double> kappas_on_spline; // 샘플링된 곡률 값 저장할 벡터
            const int num_samples_for_cost = 20; 
            double sum_abs_kappa = 0.0;
            double max_kappa_val = -numeric_limits<double>::infinity(); 
            double min_kappa_val = numeric_limits<double>::infinity();  

            for(int k = 0; k <= num_samples_for_cost; ++k){
                double t_eval = static_cast<double>(k) / num_samples_for_cost;
                // EdgeInfo에 저장된 원본 스플라인 계수 사용
                SplinePoint sp = evaluateSpline(edgeInfo.coeffs_x_orig, edgeInfo.coeffs_y_orig, t_eval, edgeInfo.spline_len, true);
                
                kappas_on_spline.push_back(sp.kappa); // 각 샘플 지점의 곡률 저장
                sum_abs_kappa += abs(sp.kappa); // 절대 곡률 합계
                max_kappa_val = max(max_kappa_val, abs(sp.kappa)); // 최대 절대 곡률
                min_kappa_val = min(min_kappa_val, abs(sp.kappa));  // 최소 절대 곡률
            }

            // 평균 곡률
            double avg_abs_kappa = sum_abs_kappa / static_cast<double>(num_samples_for_cost + 1);
            current_offline_cost += params.w_curv_avg * pow(avg_abs_kappa, 2) * edgeInfo.spline_len;

            // 피크 곡률
            double diff_peak_kappa = abs(max_kappa_val - min_kappa_val); // 최대 절대 곡률 - 최소 절대 곡률
            current_offline_cost += params.w_curv_peak * pow(diff_peak_kappa, 2) * edgeInfo.spline_len;

            // 경로 길이 비용
            current_offline_cost += params.w_length * edgeInfo.spline_len;

            // 레이싱 라인 비용 계산
            int raceline_idx_at_end_layer = -1;

            // 현재 edge의 끝 node 인덱스가 속한 레이어는 srcKey의 layer + 1
            int end_layer_idx =  (static_cast<int>(get<0>(srcKey)) + 1) % static_cast<int>(nodesPerLayer.size());
            
            // 레이싱 라인 플래그 true 인 노드 찾기
            if (end_layer_idx < nodesPerLayer.size()) {
                for (const auto& node_in_end_layer : nodesPerLayer[end_layer_idx]) {
                    if (node_in_end_layer.raceline) {
                        raceline_idx_at_end_layer = node_in_end_layer.node_idx;
                        break;
                    }
                }
            }
            
            // 레이싱 라인 노드
            if (raceline_idx_at_end_layer != -1) {
                double raceline_dist = abs(raceline_idx_at_end_layer - destIdx) * params.lat_resolution;
                current_offline_cost += min(params.w_raceline * edgeInfo.spline_len * raceline_dist,
                                                 params.w_raceline_sat * edgeInfo.spline_len);
            } else {
                // cout << "Warning: Raceline node not found for end layer " << end_layer_idx << ". Raceline cost not fully applied." << endl;
            }

            // 최종 비용 edgeInfo에 저장
            edgeInfo.offline_cost = current_offline_cost;

            processed_edges_count++;
        }
    }

}

void set_startpos(const Vector2d& pos_est, double heading_est, const Offline_Params& params, const NodeMap& nodesPerLayer,
                  Graph& graph, double max_heading_offset_rad, ITuple& out_start_key, int& out_next_idx, bool& out_of_track){
    // 레이어 0 노드 선택 > 레이어 1 노드 선택 > 엣지 존재 여부 확인 > 없으면 생성 / 있으면 유지
    
    if(nodesPerLayer.empty()) throw runtime_error("Graph not initialized.");

    // 1. 트랙 내부 여부
    const bool in_track = isPointInsideTrackBounds(pos_est.x(), pos_est.y());

    // 2. 가장 가까운 레이어 인덱스
    int l0_tmp = -1; double best = numeric_limits<double>::max();
    for(size_t i = 0; i < sampling_map[__x_ref].size(); ++i){
        double dx = pos_est.x() - sampling_map[__x_ref][i];
        double dy = pos_est.y() - sampling_map[__y_ref][i];
        double d2 = dx*dx + dy*dy;
        if(d2 < best) { best = d2; l0_tmp = (int)i; }
    }

    // l0_tmp 가 유효할 때만 헤딩 비교 가능
    bool cor_heading = false;
    if(l0_tmp >= 0 && l0_tmp < (int)sampling_map[__psi].size()){
        cor_heading = abs(normalizeAngle(heading_est - sampling_map[__psi][l0_tmp])) <= max_heading_offset_rad;
    }

    // 하나라도 실패하면 그냥 종료
    if(!in_track || !cor_heading || l0_tmp < 0 || l0_tmp >= (int)nodesPerLayer.size()){
        out_of_track = true;
        out_start_key = ITuple(-1, -1);
        out_next_idx = -1;
        return;
    }

    const int l0 = l0_tmp;

    // 3. 레이어 0 노드 인덱스
    const double ref_x = sampling_map[__x_ref][l0], ref_y = sampling_map[__y_ref][l0];
    const double nx = sampling_map[__x_normvec][l0], ny = sampling_map[__y_normvec][l0];
    const double wL = sampling_map[__width_left][l0], wR = sampling_map[__width_right][l0];
    const double alpha_v = (pos_est.x() - ref_x) * nx + (pos_est.y() - ref_y) * ny;

    const double lat = params.lat_resolution, vehW = params.veh_width;
    const int raceline_index = (int)floor((wL + sampling_map[__alpha][l0] - 0.5 * vehW) / lat);
    const double start_alpha = sampling_map[__alpha][l0] - raceline_index * lat;

    int num_nodes = (int)((wR + wL - vehW) / lat) + 1;
    num_nodes = max(num_nodes, 1);

    int node_idx0 = (int)llround((alpha_v - start_alpha) / lat);
    node_idx0 = clamp(node_idx0, 0, min(num_nodes - 1, (int)nodesPerLayer[l0].size() - 1));

    out_start_key = ITuple(l0, node_idx0);

    // 4. 레이어 1 노드
    const int l1 = (l0 + 1) % (int)nodesPerLayer.size();
    out_next_idx = clamp(node_idx0, 0, (int)nodesPerLayer[l1].size() - 1);
    // 시작 노드와 같은 횡방향 인덱스(node_idx0)를 가진 다음 레이어(l1)의 노드를 목적지로 선택
    // 트랙을 따라 자연스럽게 진행하는 경로를 가정

    // 5. 엣지 없으면 생성
    const auto& adj = graph.getAdjLists();
    const auto it = adj.find(out_start_key);
    const bool has_edge = (it != adj.end()) && (it->second.count(out_next_idx) > 0);

    if(!has_edge){
        const Node& a = nodesPerLayer[l0][node_idx0];
        const Node& b = nodesPerLayer[l1][out_next_idx];

        MatrixXd P(2, 2); P << a.x, a.y, b.x, b.y;
        VectorXd L(1); L(0) = (P.row(1) - P.row(0)).norm();

        try{
            auto res = calcSplines(P, &L, a.psi, b.psi, true);
            if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params))
                graph.addEdge(out_start_key, out_next_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
        }catch(...){}
    }
    
    out_of_track = false;
}



