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

    // 닫힌 경로 확인(np.all(np.isclose(path[0], path[-1]))와 동일)
    bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6;

    if (closed){
        // 닫힌 경로의 경우 path 마지막에 path.row(0) 추가하기 때문에 행 +1, 열 개수 변경 x
        updated_path.conservativeResize(path.rows() + 1, NoChange)
        updated_path.row(path.rows()) = path.row(0);
    }

    // 실제 스플라인이 정의될 구간의 개수 저장
    int no_splines = updated_path.rows() - 1;

    VectorXd ds;
    if (el_lengths_ptr == nullptr){
        ds = computeEuclideanDistances(updated_path); // el_lenghts_ptr 주어지지 않았다면 직접 계산
    }else{
        ds = *el_lengths_ptr // el_lengths_ptr 주어졌다면 가리키는 내용(실제 VectorXd 객체)을 ds에 복사
    }

    VectorXd scaling;
    if(use_dist_scaling){
        if(closed){
            VectorXc temp_ds(ds.size() + 1); // 현재 구간 길이 벡터보다 길이가 1 더 긴 임시 벡터 생성
            temp_ds.head(ds.size()) = ds;
            temp_ds(ds.size()) = ds(0);
            ds = temp_ds;
        }
        // ds(i) / dis(i+1) -> 현재 구간 길이 / 다음 구간 길이
        scaling = ds.head(no_splines) / ds.tail(no_splines);
    }else{ // 거리 기반 스케일링 사용 x
        scaling = VectorXd::Ones(no_splines - 1); // scaling 벡터를 모든 요소가 1인 벡터로 설정
    }
}