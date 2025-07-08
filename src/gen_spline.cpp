#include "graph_planner.hpp"
#include "config.h"

// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;          // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;          // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                 // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;// 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

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
    MatrixXd updated_path = path;

    // 닫힌 경로 확인(np.all(np.isclose(path[0], path[-1]))와 동일)
    bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6;

    if (closed){
        // 닫힌 경로의 경우 path 마지막에 path.row(0) 추가하기 때문에 행 +1, 열 개수 변경 x
        updated_path.conservativeResize(path.rows() + 1, NoChange);
        updated_path.row(path.rows()) = path.row(0);
    }

    // 실제 스플라인이 정의될 구간의 개수 저장
    int no_splines = updated_path.rows() - 1;

    VectorXd ds;
    if (el_lengths_ptr == nullptr){
        ds = computeEuclideanDistances(updated_path); // el_lenghts_ptr 주어지지 않았다면 직접 계산
    }else{
        ds = *el_lengths_ptr; // el_lengths_ptr 주어졌다면 가리키는 내용(실제 VectorXd 객체)을 ds에 복사
    }

    VectorXd scaling;
    if(use_dist_scaling){
        if(closed){
            VectorXd temp_ds(ds.size() + 1); // 현재 구간 길이 벡터보다 길이가 1 더 긴 임시 벡터 생성
            temp_ds.head(ds.size()) = ds;
            temp_ds(ds.size()) = ds(0);
            ds = temp_ds;
        }
        // ds(i) / dis(i+1) -> 현재 구간 길이 / 다음 구간 길이
        scaling = ds.head(no_splines) / ds.tail(no_splines);
    }else{ // 거리 기반 스케일링 사용 x
        scaling = VectorXd::Ones(no_splines - 1); // scaling 벡터를 모든 요소가 1인 벡터로 설정
    }

    // M 행렬 및 우변 벡터 초기화
    MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4); // 계수 행렬
    VectorXd b_x = VectorXd::Zero(no_splines * 4); // 연립 선형 방정식 M x coeffs = b에서 결과 값 담는 벡터
    VectorXd b_y = VectorXd::Zero(no_splines * 4);

    // template_M 상수 정의
    for (int i = 0;i < no_splines; ++i){
        int j = i * 4; // 현재 spline 계수의 시작 인덱스

        // 위치 제약 조건(p_i(0) = path_i, p_i(1) = path_{i+1})
        // 현재 spline 구간 i가 시작점 path_i와 끝점 path_{i+1}을 정확히 지나도록 강제
        M(j, j) = 1; // 현재 spline의 i의 a0 계수에 해당하는 열을 의미
        b_x(j) = updated_path(i, 0);
        b_y(j) = updated_path(i, 1);

        // P(1) = a_3+a_2+a_1+a_0
        M(j + 1, j) = 1; M(j + 1, j + 1) = 1, M(j + 1, j + 2) = 1; M(j + 1, j + 3) = 1;
        b_x(j + 1) = updated_path(i + 1, 0);
        b_y(j + 1) = updated_path(i + 1, 1);

        // 연속성 제약 조건(1차 미분, 2차 미분)
        if(i < no_splines - 1){ //  마지막 spline 구간을 제외한 나머지 구간에 대해 연속성 제약
            // p_i'(1) = p_{i+1}'(0) -> 1차 미분
            // p_i'(1) = a_1 + 2a_2 + 3a_3 | p_{i+1}'(0) = a_1'
            // (a_1 + 2a_2 + 3a_3)_i - (a_1)_i+1 = 0
            M(j + 2, j + 1) = 1; M(j + 2, j + 2) = 2; M(j + 2, j + 3) = 3;
            M(j + 2, j + 5) = -scaling(i); // j + 5: 다음 spline 구간 i+1의 a1 계수, scaling(i) 곱해서 미분값 스케일 맞춰줌(방정식에 들어간 - 그대로 행렬에 반영)

            // p_i''(1) = p_{i+1}''(0) -> 2차 미분
            // p_i''(1) = 2a_2 + 6a_3 | p_{i+1}''(0) = 2a_2'
            // (2a_2 + 6a_3)_i - (2a_2)_{i+1} = 0
            M(j + 3, j + 2) = 2; M(j + 3, j + 3) = 6;
            M(j + 3, j + 6) = -2 * pow(scaling(i), 2); // 제곱 -> 2번 미분에서 각 미분에 스케일링 한 번씩 적용
        }
    }

    // 경계 조건 설정(heading) | heading: 1차 미분
    if(!closed){ // 열린 경로인 경우: 시작/끝점 헤딩 고정 (psi_s, psi_e 주어짐)
        // ---Heading start point---
        // spline 시작점의 heading 조건 위해 조건 방정식 좌변 계수 설정
        // P'(t)=a_1+2a_2t+3a_3t^2 -> P'(0)=a_1
        // M의 마지막 두 행은 경계 조건을 위한 자리(마지막 spline을 제외하고 1차, 2차 미분 연속성 제약하기 때문)
        M(no_splines * 4 - 2, 1) = 1;

        // spline 시작점의 물리적 heading 계산하기 위해 첫 번째 spline 구간 길이 가져옴.
        double el_lengths_s = (el_lengths_ptr == nullptr) ? 1.0 : ds(0); // 길이 정보가 없으면 기본값으로 1.0
        
        // b_x, b_y 우변 벡터에 [시작점 헤딩(psi_s) = cos/sin 변환 값 * el_length_s] 설정
        b_x(no_splines * 4 - 2) = cos(psi_s + M_PI / 2) * el_length_s;
        b_y(no_splines * 4 - 2) = sin(psi_s + M_PI / 2) * el_length_s;


        // ---Heading end point---
        // P'(t)=a_1+2a_2t+3a_3t^2 -> P'(1)=a_1+2a_2+3a_3
        int last_spline_idx_start = 4 * (no_splines - 1);
        M(no_splines * 4 - 1, last_spline_idx_start) = 0; // a_0
        M(no_splines * 4 - 1, last_spline_idx_start + 1) = 1; // a_1 
        M(no_splines * 4 - 1, last_spline_idx_start + 2) = 2; // a_2
        M(no_splines * 4 - 1, last_spline_idx_start + 3) = 3; // a_3

        double el_lengths_e = (el_lenghts_ptr == nullptr) ? 1.0 : ds(no_splines - 1);

        // b_x, b_y 우변 벡터에 [끝점 헤딩(psi_e) = cos/sin 변환 값 * el_length_e] 설정        
        b_x(no_splines * 4 - 1) = cos(psi_e + M_PI / 2) * el_length_e;
        b_y(no_splines * 4 - 1) = sin(psi_e + M_PI / 2) * el_length_e;
    }else{ // 닫힌 경로인 경우: heading/curvature 주기 조건(첫 spline 시작 = 마지막 spline끝)
        // Heading 경계 조건
        // p_0'(0) - p_{last}'(1) = 0
        // a1_0 - a1_last - 2*a2_last - 3*a3_last = 0
        M(no_splines * 4 - 2, 1) = 1; // a1_0
        int last_spline_idx_start = 4 * (no_splines - 1);
        M(no_splines * 4 - 2, last_spline_idx_start + 1) = -scaling(no_splines - 1);
        M(no_splines * 4 - 2, last_spline_idx_start + 2) = -2 * scaling(no_splines - 1);
        M(no_splines * 4 - 2, last_spline_idx_start + 3) = -3 * scaling(no_splines - 1);
        // 이미 b_x, b_y는 0으로 설정되어 있음.

        // Curvature 경계 조건
        // p_0''(0) - p_{last}''(1) = 0
        // 2*a2_0 - 2*a2_last - 6*a3_last = 0
        M(no_splines * 4 - 1, 2) = 2; // 2a2_0
        M(no_splines * 4 - 1, last_spline_idx_start + 2) = -2 * pow(scaling(no_splines - 1), 2);
        M(no_splines * 4 - 1, last_spline_idx_start + 3) = -6 * pow(scaling(no_splines - 1), 2);
        // 이미 b_x, b_y는 0으로 설정되어 있음.
    }

    // 연립방정식 풀기(M · coeffs = b)
    // x/y_les: spline의 x, y 좌표 나타내는 모든 3차 다항식들의 계수들 나열한 1차원 벡터
    VectorXd x_les = M.colPivHouseholderQr().solve(b_x);
    VectorXd y_les = M.colPivHouseholderQr().solve(b_y);

    // 계수 행렬로 변환
    // .data(): 데이터가 저장된 메모리 시작 주소 가져옴.
    // 4: row 개수(다항식 계수 4개)
    // transpose -> 각 행이 하나의 spline 구간의 4개 계수([a0, a1, a2, a3])를 담게 됨.
    MatrixXd coeffs_x_res = Map<MatrixXd>(x_les.data(), 4, no_splines).transpose(); 
    MatrixXd coeffs_y_res = Map<MatrixXd>(y_les.data(), 4, no_splines).transpose();

    // 법선 벡터 계산
    MatrixXd normvec(no_splines, 2); // 각 spline 구간에 대한 2차원 법선 벡터
    for(int i = 0; i < no_splines; ++i){
        // i번째 스플라인 구간의 X/Y 방향 1차항 계수(a_1)
        double dx = coeffs_x_res(i, 1);
        double dy = coeffs_y_res(i, 1);
        // 계산된 접선 벡터에 대해 수직인 벡터 저장
        normvec(i, 0) = -dy;
        normvec(i, 1) = dx;
    }

    // 법선 벡터 정규화
    VectorXd norms = normvec.rowwise().norm();
    MatrixXd normvec_normalized_res(no_splines, 2); // normvec 행렬의 각 행에 대해 유클리드 노름 계산, norm 벡터에 저장
    for(int i = 0; i < no_splines; ++i){
        if(norms(i) > 1e-9){ // 0으로 나누는 오류 방지
            normvec_normalized_res.row(i) = normvec.row(i) / norms(i); // 단위 벡터 생성
        }else{
            normvec_normalized_res.row(i).setZero(); // 0으로 나누는 경우 0 벡터로 설정
        }        
    }

    return {coeffs_x_res, coeffs_y_res, M, normvec_normalized_res, ds}; // 각 구간의 실제 길이인 ds 값도 함께 반환
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

    // lateral_offset이 트랙의 유효한 횡방향 범위 내에 있는지 직접적으로 확인하는 최종 단계
    if(lateral_offset >= -width_left && lateral_offset <= width_right){
        return true;
    }else{
        return false;
    }
}

// 하나의 스플라인 구간이 유효한 경로로 사용될 수 있는지 검사
bool checkSplineValidity(const RowVector4d coeff_x, const RowVector4d& coeff_y, double ds_current,
                         const Offline_Params& params){ // const MatrixXd& track_bounds
    // spline 경로 샘플링
    const int num_samples = 20;

    for(int k = 0; k <= num_samples; ++k){
        double t_eval = static_cast<double>(k) / num_samples; // 파라미터 t 값을 균등하게 분할
        SplinePoint sp = evaluateSpline(coeff_x, coeff_y, t_eval, ds_current, true);

        // 트랙 경계 벗어나는지 확인
        if(!isPointInsideTrackBounds(sp.x, sp.y)){
            return false;
        }

        // 곡률 제약 조건 확인
        // std::abs(sp.kappa): 현재 spline 점에서의 kappa 절댓값
        // 1.0 / params.VEH_TURN: 차량이 허용하는 최대 곡률 (최소 회전 반경의 역수)
        if(std::abs(sp.kappa) > (1.0 / params.VEH_TURN)){
            return false;
        }
    }

    return true;
}

void generateGraphEdges(Graph& graph, const NodeMap& nodesPerLayer, const Offline_Params& params){
    const size_t num_layers = nodesPerLayer.size(); // 총 layer의 개수

    // current_layer_idx에서 next_layer_idx로의 연결 시도
    for(size_t current_layer_idx = 0; current_layer_idx < num_layers; ++current_layer_idx){
        // 다음 layer의 인덱스 계산
        size_t next_layer_idx = (current_layer_idx + 1) % num_layers; // 닫힌 트랙에서 마지막 layer 이후 첫 번째 layer로 돌아가게 함.

        // 현재 layer와 다음 layer의 node 가져옴.
        const auto& current_nodes_in_layer = nodesPerLayer[current_layer_idx];
        const auto& next_nodes_in_layer = nodesPerLayer[next_layer_idx];

        // 현재 layer의 모든 current_node에 대해 수행
        for(const auto& current_node : current_nodes_in_layer){
            for(const auto& next_node : next_nodes_in_layer){
                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y, next_node.x, next_node.y; // 첫 번째 행: 현재 node의 (x, y) | 두 번째 행: 다음 node의 (x, y)

                // spline 시작/끝 heading으로 각 node의 psi값 사용
                double psi_s = current_node.psi;
                double psi_e = next_node.psi;

                // spline 구간의 유클리드 길이
                VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

                SplineResult res;
                try{
                    // spline_path의 두 점을 연결하는 spline 계수 계산
                    res = calcSplines(spline_path, &el_lengths, psi_s, psi_e, true);
                }catch(const std::exception& e){
                    continue; // 유효하지 않을 경우(예외 발생) 건너뜀.
                }

                // spline 유효성 검사
                if(checkSplineValidity(res.coeffs_x.row(0), res.coeffs_y.row(0), res.ds(0), params, MatrixXd())){
                    ITuple src_key(current_node.layer_idx, current_node.node_idx);
                    graph.addEdge(src_key, next_node.node_idx); // 검사 통과 시 그래프에 edge 추가(src_key: 특정 node를 고유하게 식별하는 key 역할, tuple)
                }
            }
        }
    }
}

int main(){

}