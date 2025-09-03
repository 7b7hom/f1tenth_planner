#include "graph_planner.hpp"
#include "config_modena.h"
#include <fstream>
#include <limits>
#include <utility>


// 앞뒤 공백/CR/LF 제거 + UTF-8 BOM 제거
static inline std::string trim_and_debom(std::string s) {
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };

    // ltrim
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    // rtrim (공백/탭/CR/LF 모두 제거)
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());

    // UTF-8 BOM 제거
    if (s.size() >= 3 &&
        (unsigned char)s[0] == 0xEF &&
        (unsigned char)s[1] == 0xBB &&
        (unsigned char)s[2] == 0xBF) {
        s.erase(0,3);
    }
    return s;
}

/// ----------------------------
//  I/O
// ----------------------------

// CSV를 읽어서 DMap으로 변경
/*DMap readDMapFromCSV(const string& pathname) {
    DMap map;
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();
    for (const auto& label : labels) map[label] = csv.GetColumn<double>(label);

    return map;
}*/

DMap readDMapFromCSV(const string& pathname) {
    DMap map;
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& raw : labels) {
        string clean = trim_and_debom(raw);         // ← 정규화된 키
        map[clean] = csv.GetColumn<double>(raw);    // 값은 원래 라벨로 읽어도 됨
    }

    // 필요시 한 번만 켜서 실제 바이트 확인
    // dump_available_keys_verbose(map);

    return map;
}


// DMap을 CSV에 작성
size_t writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',') {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    size_t num_cols = map.size();
    size_t num_rows = map.empty() ? 0 : map.begin()->second.size();

    // header
    size_t i = 0;
    for (const auto& [key, _] : map) {
        file << key; if (++i != num_cols) file << delimiter;
    }
    file << '\n';

    // rows
    for (size_t row = 0; row < num_rows; ++row) {
        size_t j = 0;
        for (const auto& [_, col] : map) {
            file << col[row]; if (++j != num_cols) file << delimiter;
        }
        file << '\n';
    }

    return num_rows;
}

// Debug용 함수: map의 columns, rows 개수 print
void map_size(DMap& map) {
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

XYPair makeXYFromOffset(const DVector& x_ref, const DVector& y_ref, const DVector& nx, const DVector& ny, const DVector& offset){
    const size_t n = x_ref.size();
    if(y_ref.size()!=n || nx.size()!=n || ny.size()!=n || offset.size()!=n){
        throw runtime_error("makeXYFromOffset: input size mismatch");
    }
    DVector xs(n), ys(n);
    for(size_t i = 0; i < n; ++i){
        xs[i] = x_ref[i] + nx[i] * offset[i];
        ys[i] = y_ref[i] + ny[i] * offset[i];
    }
    return { move(xs), move(ys) };
}

XYPair computeBoundRight(const DMap& src){
    return makeXYFromOffset(src.at(__x_ref), src.at(__y_ref), 
                            src.at(__x_normvec), src.at(__y_normvec), 
                            src.at(__width_right));
}

XYPair computeBoundLeft(const DMap& src){
    const auto& wL = src.at(__width_left);
    DVector neg_wL(wL.size());
    for(size_t i = 0; i < wL.size(); ++i) neg_wL[i] = -wL[i];

    return makeXYFromOffset(src.at(__x_ref), src.at(__y_ref), 
                            src.at(__x_normvec), src.at(__y_normvec), 
                            neg_wL);
}

XYPair computeBoundRace(const DMap& src){
    return makeXYFromOffset(src.at(__x_ref), src.at(__y_ref),
                            src.at(__x_normvec), src.at(__y_normvec), 
                            src.at(__alpha));
}

DVector computeDeltaS(const DVector& s, bool closed){
    const size_t n = s.size();
    DVector ds(n, 0.0);
    for(size_t i = 0; i + 1 < n; ++i) ds[i] = s[i+1] - s[i];

    if(closed && n >= 2) ds[n-1] = 0.0;
    return ds;
}

// 레이싱 라인에서 경로 계획을 위한 layer가 될 지점들을 샘플링
IVector samplePointsFromRaceline(const DVector& kappa, const DVector& dist,
                              double d_curve, double d_straight, double curve_th) {
    IVector idx;
    const size_t n = kappa.size();
    double cur = 0.0, next = 0.0, next_min = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if ((cur + dist[i]) > next_min && fabs(kappa[i]) > curve_th) next = cur;
        if ((cur + dist[i]) > next) {
            idx.push_back((int)i);
            next += (fabs(kappa[i]) < curve_th) ? d_straight : d_curve;
            next_min = cur + d_curve;
        }
        cur += dist[i];
    }
    return idx;
}

// 주어진 각도를 -PI에서 PI 사이의 값으로 정규화
double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// 주어진 X, Y 좌표 벡터를 기반으로 각 지점에서의 헤딩(Heading, 진행 방향 각도)을 계산
DVector calcHeading(const DVector& x, const DVector& y) {
    size_t N = x.size();
    DVector psi(N);
    for (size_t i = 0; i < N; ++i) {
        double dx, dy;
        if(i != N-1) { dx = x[i+1] - x[i]; dy = y[i+1] - y[i]; }
        else { dx = x[0] - x[i], dy = y[0] - y[i]; }
        psi[i] = normalizeAngle(atan2(dy, dx) - M_PI_2);
    }
    return psi;
}

#if 0
// 샘플링된 레이어마다 경로 계획을 위한 Node(차량이 횡방향으로 이동 가능한 위치들) 생성
NodeMap genNode(const DMap& sampled_map, const double veh_width, float lat_resolution) {
    const size_t N = sampled_map.at(__alpha).size();
    NodeMap nodes(N);

    for (size_t i = 0; i < N; ++i){ // 각 레이어(층)에 대해 반복
        const double wL = sampled_map.at(__width_left)[i];
        const double wR = sampled_map.at(__width_right)[i];
        const double alpha = sampled_map.at(__alpha)[i];

        const int raceline_index = (int)floor((wL + alpha - veh_width / 2.0) / lat_resolution);
        const double start_alpha = alpha - raceline_index * lat_resolution;
        const int num_nodes = (int)((wL + wR - veh_width) / lat_resolution) + 1;

        nodes[i].resize(max(num_nodes, 1));
        
        Vector2d ref_xy(sampled_map.at(__x_ref)[i], sampled_map.at(__y_ref)[i]);
        Vector2d norm_vec(sampled_map.at(__x_normvec)[i], sampled_map.at(__y_normvec)[i]);
        
        for (int idx = 0; idx < (int)nodes[i].size(); ++idx) { 
            double alphaI = start_alpha + idx * lat_resolution; // 현재 노드의 횡방향 오프셋 계산
            Vector2d node_pos = ref_xy + alphaI * norm_vec; // 노드의 (x, y) 좌표 계산 (여기서 선언)

            Node node; // 각 노드 인스턴스를 이 루프 안에서 새로 생성하여 초기화 문제를 방지
            
            // 필수 멤버 초기화 및 할당
            node.layer_idx = i; // 현재 레이어 인덱스
            node.node_idx = idx; // 노드 인덱스
            node.x = node_pos.x();
            node.y = node_pos.y();
            node.kappa = sampled_map.at(__kappa)[i];
            node.raceline = (idx == raceline_index);

            // --- 노드의 헤딩(psi) 계산 (보간) ---
            double psi_interp;
            if (idx < raceline_index) { 
                double bl = sampled_map.at(__psi_bound_l)[i];
                double pr = sampled_map.at(__psi)[i];
                if(abs(bl - pr) >= M_PI){
                    bl += 2*M_PI * (bl < 0);
                    pr += 2*M_PI * (pr < 0);
                }
                psi_interp = bl + (pr - bl) * (double)(idx + 1) / max(raceline_index, 1);
            } else if (idx == raceline_index) { 
                psi_interp = sampled_map.at(__psi)[i];
            }
            else { 
                int remain = max((int)(nodes[i].size()) - raceline_index - 1, 1);
                double t = double(idx - raceline_index) / remain;
                psi_interp = sampled_map.at(__psi)[i] + t * (sampled_map.at(__psi_bound_r)[i] - sampled_map.at(__psi)[i]);
            }
            node.psi = normalizeAngle(psi_interp);
            
            nodes[i][idx] = node; // <-- 생성된 노드 인스턴스를 NodeMap에 할당
        }
    }
    return nodes;
}
#endif

// --- genNode 3단계 체인 ---

// 1. layer 파라미터 전처리
vector<LayerParams> computeLayerParams(const DMap& map, double vehW, float lat){
    const size_t N = map.at(__alpha).size();
    vector<LayerParams> result; result.reserve(N);

    for(size_t i = 0; i < N; ++i){
        LayerParams p;
        p.layer_idx = (int)i;
        p.wL = map.at(__width_left)[i];
        p.wR = map.at(__width_right)[i];
        p.alpha = map.at(__alpha)[i];
        p.ref_xy = { map.at(__x_ref)[i], map.at(__y_ref)[i] };
        p.norm_vec = { map.at(__x_normvec)[i], map.at(__y_normvec)[i] };

        p.raceline_index = (int)floor((p.wL + p.alpha - 0.5*vehW) / lat);
        p.start_alpha = p.alpha - p.raceline_index * lat;
        p.num_nodes = max((int)((p.wL + p.wR - vehW) / lat) + 1, 1);

        result.push_back(move(p));
    }
    return result;
}

// 2. 좌표/기본 필드로 Node 그리드 생성
NodeMap buildNodeGrid(const vector<LayerParams>& layers, const DMap& map, float lat){
    NodeMap nodes(layers.size());

    for (size_t i = 0; i < layers.size(); ++i){
        const auto& L = layers[i];
        nodes[i].resize(L.num_nodes);

        for (int idx = 0; idx < L.num_nodes; ++idx){
            double alphaI = L.start_alpha + idx * lat;
            Vector2d pos = L.ref_xy + alphaI * L.norm_vec;

            Node n;
            n.layer_idx = (int)i;
            n.node_idx = idx;
            n.x = pos.x();
            n.y = pos.y();
            n.kappa = map.at(__kappa)[i];
            n.raceline = (idx == L.raceline_index);

            nodes[i][idx] = n;
        }
    }
    return nodes;
}

// 3. 노드 헤딩 보간, 채우기
void fillNodeHeadings(NodeMap& nodes, const vector<LayerParams>& layers, const DMap& map){
    for (size_t i = 0; i < nodes.size(); ++i){
        if (nodes[i].empty()) continue;

        const double psi_race = map.at(__psi)[i];
        const double psi_left = map.at(__psi_bound_l)[i];
        const double psi_right = map.at(__psi_bound_r)[i];

        int R = layers[i].raceline_index;
        int N = static_cast<int>(nodes[i].size());

        R = clamp(R, 0, max(N - 1, 0)); // 0 <= R < N 범위

        for (int idx = 0; idx < N; ++idx){
            double psi_interp;

            if (idx < R){
                // 왼쪽 경계 > 레이스 라인 보간
                const int denom = max(R, 1);
                const double t = static_cast<double>(idx + 1) / denom;
                const double diff = normalizeAngle(psi_race - psi_left);
                psi_interp = normalizeAngle(psi_left + t * diff);
            } else if (idx == R){
                psi_interp = psi_race;
            } else {
                // 레이스 라인 > 오른쪽 경계 보간
                const int remain = max(N - R - 1, 1);
                const double t = static_cast<double>(idx - R) / remain;
                const double diff = normalizeAngle(psi_right - psi_race);
                psi_interp = normalizeAngle(psi_race + t * diff);
            }

            nodes[i][idx].psi = psi_interp;
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

// 주어진 (x, y) 점이 트랙 경계 내부에 있는지(왼, 오 경계 사이에 위치하는지) 확인하는 함수 (전역 sampled_map 사용)
bool isPointInsideTrackBounds(const DMap& sampled_map, double x, double y){
    if(sampled_map.empty() || sampled_map.begin()->second.empty()) return false;

    // (x, y) 점에서 가장 가까운 sampled_map의 기준선 인덱스 찾기
    double best = numeric_limits<double>::max(); // 가능한 가장 큰 수로 초기값 설정
    int idx = -1; // 유효하지 않은 값으로 초기값 설정
    
    const auto& xr = sampled_map.at(__x_ref);
    const auto& yr = sampled_map.at(__y_ref);
    for(size_t i = 0; i < xr.size(); ++i){
        double dx = x - xr[i], dy = y - yr[i];
        double d2 = dx*dx + dy*dy;
        if(d2 < best) { best = d2, idx = (int)i; }
    }
    if(idx < 0) return false;

    // 가장 가까운 기준선 인덱스의 정보 가져오기
    const double ref_x = sampled_map.at(__x_ref)[idx]; // __x_ref, __y_ref: 기준선 좌표
    const double ref_y = sampled_map.at(__y_ref)[idx];
    const double norm_x = sampled_map.at(__x_normvec)[idx]; // __x_normvec, __y_normvec: 기준선 법선 벡터
    const double norm_y = sampled_map.at(__y_normvec)[idx];
    const double width_L = sampled_map.at(__width_left)[idx]; // __width_left, __width_right: 기준선으로부터의 좌우 폭
    const double width_R = sampled_map.at(__width_right)[idx];

    // spline 점의 기준선 법선 방향 횡방향 오프셋 계산
    // -> (x, y) 점이 기준선으로부터 법선 벡터 방향으로 얼마나 떨어져 있는가
    const double lateral_offset = (x - ref_x) * norm_x + (y - ref_y) * norm_y;
    
    // lateral_offset이 트랙의 유효한 횡방향 범위 내에 있는지 직접적으로 확인하는 최종 단계
    const double epsilon = 1e-6;
    return (lateral_offset >= -width_L - epsilon) && (lateral_offset <= width_R + epsilon);
}

// 하나의 스플라인 구간이 유효한 경로로 사용될 수 있는지 검사
bool checkSplineValidity(const RowVector4d& coeff_x, const RowVector4d& coeff_y, double ds_current,
                         const Offline_Params& params, const DMap& sampled_map){ // const MatrixXd& track_bounds
    // spline 경로 샘플링
    const int num_samples = 10;

    double max_allowed_kappa = 30.0 / params.veh_turn;

    for(int k = 0; k <= num_samples; ++k){
        double t_eval = static_cast<double>(k) / num_samples; // 파라미터 t 값을 균등하게 분할
        SplinePoint sp = evaluateSpline(coeff_x, coeff_y, t_eval, ds_current, true);

        // 트랙 경계 벗어나는지 확인
        if(!isPointInsideTrackBounds(sampled_map, sp.x, sp.y)){
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

#if 0
Graph genEdges(const NodeMap& nodesPerLayer, const Offline_Params& params, const DMap& sampled_map){
    Graph graph(true);
    const size_t layers = nodesPerLayer.size();

    auto& adj = graph.getAdjLists_mutable();
    for (size_t l = 0; l < layers; ++l) {
        for (size_t j = 0; j < nodesPerLayer[l].size(); ++j) {
            adj[ ITuple((int)l, (int)j) ]; // touch: 빈 map 생성
        }
    }

    // raceline 추종 엣지 살리기
    for(size_t currLayer = 0; currLayer < layers; ++currLayer){
        size_t nextLayer = (currLayer + 1) % layers;

        // 현재 layer의 레이스 라인 노드 찾기
        const Node* currN = nullptr;
        for(const auto& nd : nodesPerLayer[currLayer]) if(nd.raceline) { currN = &nd; break; }

        // 다음 layer의 레이스 라인 노드 찾기
        const Node* nextN = nullptr;
        for(const auto& nd : nodesPerLayer[nextLayer]) if(nd.raceline) { nextN = &nd; break; }
        
        if(!currN || !nextN) continue;

        Eigen::MatrixXd P(2,2); P << currN->x, currN->y, nextN->x, nextN->y;
        Eigen::VectorXd lengths(1); lengths(0) = (P.row(1) - P.row(0)).norm();

        SplineResult res = calcSplines(P, &lengths, currN->psi, nextN->psi, true);
        if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, sampled_map)){
            graph.addEdge(ITuple(currN->layer_idx, currN->node_idx), nextN->node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
        }
    }

    // layer 순회
    for(size_t currLayer = 0; currLayer < layers; ++currLayer){
        size_t nextLayer = (currLayer + 1) % layers;

        const auto& currNodes = nodesPerLayer[currLayer];
        const auto& nextNodes = nodesPerLayer[nextLayer];

        // 현재 layer의 node 순회
        for(const auto& currNode : currNodes){
            // lat_steps 로직
            int refDestIdx = clamp(currNode.node_idx, 0, (int)nextNodes.size()-1);
            const Node& refEndNode = nextNodes[refDestIdx];

            // 현재 layer와 다음 layer 기준 노드 사이의 거리 계산
            double dist = (Vector2d(refEndNode.x, refEndNode.y) - Vector2d(currNode.x, currNode.y)).norm();

            double ratio = (params.curve_thr > 1e-9) ? min(abs(currNode.kappa)/params.curve_thr, 2.0) : 0.0;
            
            // factor -> 직선보다 곡선에서 더 넓은 노드 탐색을 가능하게 하는 계수
            double factor = 1.0 / (1.0 + 0.5 * ratio); // 곡률 높을수록 lat_steps 줄이기

            // 곡률이 높고 거리가 멀수록 더 많은 노드를 살펴봄
            int lat_steps = static_cast<int>(round(factor * dist * params.lat_offset / params.lat_resolution));
            lat_steps = min(lat_steps, (int)params.max_lat_steps);

            for(int destNode = max(0, refDestIdx - lat_steps); destNode <= min(static_cast<int>(nextNodes.size() - 1), refDestIdx + lat_steps); ++destNode){
                const Node& nextNode = nextNodes[destNode];

                MatrixXd P(2, 2); P << currNode.x, currNode.y, nextNode.x, nextNode.y;
                VectorXd lengths(1); lengths(0) = (P.row(1) - P.row(0)).norm();

                SplineResult res = calcSplines(P, &lengths, currNode.psi, nextNode.psi, true);
                
                if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, sampled_map)){
                    ITuple src_key(currNode.layer_idx, currNode.node_idx);
                    graph.addEdge(src_key, nextNode.node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0)); 
                    //cout << "SPLINE PASSED!!!! from (" << current_node.layer_idx << "," << current_node.node_idx
                         //<< ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }else{
                    //cout << "SPLINE REJECTED from (" << current_node.layer_idx << "," << current_node.node_idx
                         //<< ") to (" << next_layer_idx << "," << next_node.node_idx << ")" << "\n" << endl;
                }
            }
        }
    }
    return graph;
}
#endif

// --- genEdge 3단계 체인 ---

// 1. 빈 그래프/노드 키 초기화
Graph makeEmptyGraph(const NodeMap& nodesPerLayer){
    Graph graph(true);

    auto& adj = graph.getAdjLists_mutable();
    for (size_t l = 0; l < nodesPerLayer.size(); ++l){
        for (size_t j = 0; j < nodesPerLayer[l].size(); ++j) adj[ ITuple((int)l, (int)j) ];
    }
    return graph;
}

// 2. 레이스라인 edge 추가
void addRacelineEdges(Graph& graph, const NodeMap& nodes, const Offline_Params& params, const DMap& map){
    const size_t L = nodes.size();
    for(size_t curr = 0; curr < L; ++curr){
        size_t next = (curr + 1) % L;

        const Node* currN = nullptr; for(const auto& nd : nodes[curr]) if(nd.raceline) { currN = &nd; break; }
        const Node* nextN = nullptr; for(const auto& nd : nodes[next]) if(nd.raceline) { nextN = &nd; break; }
        if(!currN || !nextN) continue;

        MatrixXd P(2,2); P << currN->x, currN->y, nextN->x, nextN->y;
        VectorXd len(1); len(0) = (P.row(1) - P.row(0)).norm();

        auto res = calcSplines(P, &len, currN->psi, nextN->psi, true);
        if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, map)){
            graph.addEdge(ITuple(currN->layer_idx, currN->node_idx), nextN->node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
        }
    }
}

// 3. 일반 후보 엣지 추가
void addCandidateEdges(Graph& graph, const NodeMap& nodes, const Offline_Params& params, const DMap& map){
    const size_t L = nodes.size();
    for(size_t curr = 0; curr < L; ++curr){
        size_t next = (curr + 1) % L;

        const auto& currNodes = nodes[curr];
        const auto& nextNodes = nodes[next];

        for(const auto& currN : currNodes){
            int ref = clamp(currN.node_idx, 0, (int)nextNodes.size() - 1);
            const Node& reffN = nextNodes[ref];

            double dist = (Vector2d(reffN.x, reffN.y) - Vector2d(currN.x, currN.y)).norm();
            double ratio = (params.curve_thr > 1e-9) ? min(abs(currN.kappa)/params.curve_thr, 2.0) : 0.0;
            double factor = 1.0 / (1.0 + 0.5*ratio);
            int lat_steps = (int)round(factor * dist * params.lat_offset / params.lat_resolution);
            lat_steps = min(lat_steps, (int)params.max_lat_steps);

            for(int idx = max(0, ref - lat_steps); idx <= min((int)nextNodes.size() - 1, ref + lat_steps); ++idx){
                const Node& nextN = nextNodes[idx];
                
                MatrixXd P(2,2); P << currN.x, currN.y, nextN.x, nextN.y;
                VectorXd len(1); len(0) = (P.row(1) - P.row(0)).norm();

                auto res = calcSplines(P, &len, currN.psi, nextN.psi, true);
                if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, map)){
                    graph.addEdge(ITuple(currN.layer_idx, currN.node_idx), nextN.node_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
                }
            }
        }
    }
}

Graph prune_graph(Graph graph, int num_layers, bool closed) {
    int totalRemoved = 0;
    int iter = 0;

    while (true) {
        int removed = 0;
        set<pair<ITuple, int>> toRemove;
        vector<ITuple> keys;

        for (const auto& [key, _] : graph.getAdjLists()) keys.push_back(key);
        
        for (const auto& key : keys) {
            int layer = get<0>(key);
            int node = get<1>(key);

            if (!closed && (layer == 0 || layer == num_layers - 1)) continue;

            IVector children;
            vector<ITuple> parents;

            try {
                graph.getChildIdx(key, children);
            } catch (...) {
                children.clear();
            }

            graph.getParentNode(layer, node, parents, num_layers);

            if (parents.empty()) { for (int child : children) toRemove.insert({key, child}); }
            else if (children.empty()) { for (const auto& parent : parents) toRemove.insert({parent, node}); }
        }

        for (const auto& edge : toRemove) { graph.removeEdge(edge.first, edge.second); ++removed; }
                cout << "Iteration " << iter << ": Removed " << removed << " edges." << endl;

        if (removed == 0) break;

        totalRemoved += removed;
        iter++;
    }
    cout << "Total edges removed during pruning: " << totalRemoved << endl;
    return graph;
}

Graph gen_offline_cost(Graph graph, const Offline_Params& params, const NodeMap& nodesPerLayer){
    for(auto& [srcKey, destMap] : graph.getAdjLists_mutable()){ // mutable 접근자 통해 EdgeInfo 수정
        for(auto& [destIdx, edgeInfo] : destMap){ // EdgeInfo를 참조로 받아 직접 수정
            double cost = 0.0;
            const int N = 20;
            double kappaSum = 0.0;
            double kMax = -numeric_limits<double>::infinity(); 
            double kMin = numeric_limits<double>::infinity();

            for(int k = 0; k <= N; ++k){
                double t= static_cast<double>(k) / N;
                // EdgeInfo에 저장된 원본 스플라인 계수 사용
                SplinePoint sp = evaluateSpline(edgeInfo.coeffs_x_orig, edgeInfo.coeffs_y_orig, t, edgeInfo.spline_len, true);
                
                kappaSum += abs(sp.kappa); // 절대 곡률 합계
                kMax = max(kMax, abs(sp.kappa)); // 최대 절대 곡률
                kMin = min(kMin, abs(sp.kappa));  // 최소 절대 곡률
            }

            // 평균 곡률
            double kappaAvg = kappaSum / static_cast<double>(N + 1);
            cost += params.w_curv_avg * pow(kappaAvg, 2) * edgeInfo.spline_len;

            // 피크 곡률
            cost += params.w_curv_peak * pow((kMax - kMin), 2) * edgeInfo.spline_len;

            // 경로 길이 비용
            cost += params.w_length * edgeInfo.spline_len;

            // -- 레이싱 라인 비용 계산 --
            int racelineIdx = -1;

            // 현재 edge의 끝 node 인덱스가 속한 레이어는 srcKey의 layer + 1
            int endLayer =  (static_cast<int>(get<0>(srcKey)) + 1) % static_cast<int>(nodesPerLayer.size());
            for (const auto& node : nodesPerLayer[endLayer]) if (node.raceline) { racelineIdx = node.node_idx; break; }
            if (racelineIdx != -1) {
                double dist = abs(racelineIdx - destIdx) * params.lat_resolution;
                cost += min(params.w_raceline * edgeInfo.spline_len * dist, params.w_raceline_sat * edgeInfo.spline_len);
            }
            edgeInfo.offline_cost = cost;
        }
    }
    return graph;
}

StartPosResult set_startpos(const Vector2d& pos_est, double heading_est, const Offline_Params& params, const NodeMap& nodesPerLayer,
                  const DMap& sampled_map, double max_heading_offset_rad){
    // 레이어 0 노드 선택 > 레이어 1 노드 선택 > 엣지 존재 여부 확인 > 없으면 생성 / 있으면 유지
    StartPosResult out{};
    out.out_of_track = true;
    out.start_key = ITuple(-1, -1);
    out.next_idx = -1;

    if (nodesPerLayer.empty()) return out;

    // 1. 트랙 내부 여부
    const bool in_track = isPointInsideTrackBounds(sampled_map, pos_est.x(), pos_est.y());

    // 2. 가장 가까운 레이어 인덱스
    int l0tmp = -1; double best = numeric_limits<double>::max();
    for(size_t i = 0; i < sampled_map.at(__x_ref).size(); ++i){
        double dx = pos_est.x() - sampled_map.at(__x_ref)[i];
        double dy = pos_est.y() - sampled_map.at(__y_ref)[i];
        double d2 = dx*dx + dy*dy;
        if(d2 < best) { best = d2; l0tmp = (int)i; }
    }

    // l0tmp 가 유효할 때만 헤딩 비교 가능
    bool corHeading = false;
    if(l0tmp >= 0 && l0tmp < (int)sampled_map.at(__psi).size()){
        corHeading = abs(normalizeAngle(heading_est - sampled_map.at(__psi)[l0tmp])) <= max_heading_offset_rad;
    }

    // 하나라도 실패하면 그냥 종료
    if(!in_track || !corHeading || l0tmp < 0 || l0tmp >= (int)nodesPerLayer.size()) return out;

    const int l0 = l0tmp;

    // 3. 레이어 0 노드 인덱스
    const double ref_x = sampled_map.at(__x_ref)[l0], ref_y = sampled_map.at(__y_ref)[l0];
    const double nx = sampled_map.at(__x_normvec)[l0], ny = sampled_map.at(__y_normvec)[l0];
    const double wL = sampled_map.at(__width_left)[l0], wR = sampled_map.at(__width_right)[l0];
    const double alpha_v = (pos_est.x() - ref_x) * nx + (pos_est.y() - ref_y) * ny;

    const double lat = params.lat_resolution, vehW = params.veh_width;
    const int racelineIdx = (int)floor((wL + sampled_map.at(__alpha)[l0] - 0.5 * vehW) / lat);
    const double startAlpha = sampled_map.at(__alpha)[l0] - racelineIdx * lat;

    int numNodes = (int)((wR + wL - vehW) / lat) + 1;
    numNodes = max(numNodes, 1);

    int Idx0 = (int)llround((alpha_v - startAlpha) / lat);
    Idx0 = clamp(Idx0, 0, min(numNodes - 1, (int)nodesPerLayer[l0].size() - 1));

    out.start_key = ITuple(l0, Idx0);

    // 4. 레이어 1 노드
    const int l1 = (l0 + 1) % (int)nodesPerLayer.size();
    out.next_idx = clamp(Idx0, 0, (int)nodesPerLayer[l1].size() - 1);
    // 시작 노드와 같은 횡방향 인덱스(node_idx0)를 가진 다음 레이어(l1)의 노드를 목적지로 선택
    // 트랙을 따라 자연스럽게 진행하는 경로를 가정

    out.out_of_track = false;

    return out;
}

Graph ensure_edge_between(Graph graph, const NodeMap& nodesPerLayer, const Offline_Params& params,
                          const ITuple& start_key, int next_idx, const DMap& sampled_map){
    const auto& adj = graph.getAdjLists();
    bool has = false;
    if (auto it = adj.find(start_key); it != adj.end()){
        has = (it->second.count(next_idx) > 0);
    }
    if (has) return graph;

    int l0 = get<0>(start_key);
    int i0 = get<1>(start_key);
    int l1 = (l0 + 1) % (int)nodesPerLayer.size();

    const Node& n0 = nodesPerLayer[l0][i0];
    const Node& n1 = nodesPerLayer[l1][next_idx];

    MatrixXd P(2,2); P << n0.x, n0.y, n1.x, n1.y;
    VectorXd L(1); L(0) = (P.row(1) - P.row(0)).norm();

    auto res = calcSplines(P, &L, n0.psi, n1.psi, true);
    if (checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, sampled_map)) {
        graph.addEdge(start_key, next_idx, res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0));
    }

    return graph;
}

void dump_available_keys(const DMap& m) {
    cerr << "[DMap keys] ";
    for (const auto& [k, _] : m) cerr << k << " ";
    cerr << endl;
}

void assert_has_keys(const DMap& m, const vector<string>& keys, const char* where) {
    for (const auto& k : keys) {
        if (m.find(k) == m.end()) {
            cerr << "\n[ERROR] Missing key '" << k << "' in " << where << endl;
            dump_available_keys(m);
            throw runtime_error(string("Missing key: ") + k);
        }
    }
}
