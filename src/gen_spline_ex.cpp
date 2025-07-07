#include "graph_planner.hpp"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;          // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;          // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                 // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;// 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

// 경로 행렬의 각 점 사이의 유클리드 거리 계산 함수
VectorXd computeEuclideanDistances(const MatrixXd& path) {
    int N = path.rows() - 1;    // 구간 개수는 점 개수 - 1
    VectorXd dists(N);          // 거리 저장 벡터
    for (int i = 0; i < N; ++i) {
        dists(i) = (path.row(i + 1) - path.row(i)).norm();  // 두 점 사이 거리 계산 (유클리드 노름)
    }
    return dists;
}

// 스플라인 계수를 계산하는 함수
// path: 입력 경로 (Nx2 행렬)
// el_lengths: 각 구간 길이 벡터 (옵션)
// psi_s, psi_e: 시작과 끝 방향 (라디안), NAN이면 무시
// use_dist_scaling: 구간 길이 스케일링 사용 여부 (현재 코드에서는 미사용)
SplineResult calcSplines(const MatrixXd& path, const VectorXd* el_lengths = nullptr,
                         double psi_s = NAN, double psi_e = NAN,
                         bool use_dist_scaling = true) {
    MatrixXd updated_path = path; 
    // 경로가 닫혀있는지 확인 (첫 점과 마지막 점이 거의 같으면 닫힘)
    bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6;

    if (closed) {
        // 닫힌 경로라면 시작점과 동일한 점을 끝에 한 번 더 추가하여 폐곡선 유지
        updated_path.conservativeResize(path.rows() + 1, NoChange);
        updated_path.row(path.rows()) = path.row(1);
    }

    int no_splines = updated_path.rows() - 1;  // 스플라인 구간 개수

    // 각 구간 길이 벡터 ds 계산 또는 입력값 사용
    VectorXd ds;
    if (el_lengths == nullptr) {
        ds = computeEuclideanDistances(updated_path);  // 기본 거리 계산 함수 호출
    } else {
        ds = *el_lengths;  // 외부에서 받은 길이 벡터 사용
    }

    // 시스템 행렬과 벡터 초기화
    MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
    VectorXd b_x = VectorXd::Zero(no_splines * 4);
    VectorXd b_y = VectorXd::Zero(no_splines * 4);

    // 각 스플라인 구간에 대해 위치 제약 조건 (시작점, 끝점 위치 맞춤)
    for (int j = 0; j < no_splines; ++j) {
        int row0 = j * 4;

        // 시작점 위치: p_j(0) = path[j]
        M(row0, row0) = 1;
        b_x(row0) = updated_path(j, 0);
        b_y(row0) = updated_path(j, 1);

        // 구간 길이 d에 대한 다항식 벡터 T = [1, d, d^2, d^3]
        RowVector4d T;
        double d = ds(j);
        T << 1, d, pow(d, 2), pow(d, 3);

        // 끝점 위치: p_j(d) = path[j+1]
        M(row0 + 1, row0) = T(0);
        M(row0 + 1, row0 + 1) = T(1);
        M(row0 + 1, row0 + 2) = T(2);
        M(row0 + 1, row0 + 3) = T(3);

        b_x(row0 + 1) = updated_path(j + 1, 0);
        b_y(row0 + 1) = updated_path(j + 1, 1);
    }

    // 1차 미분 연속성 제약: 인접 스플라인 접점에서 1차 미분값이 같아야 함
    for (int j = 0; j < no_splines - 1; ++j) {
        int row = 2 * no_splines + j;  // 행 번호
        double d = ds(j);

        // 1차 미분 다항식 벡터: p'(d) = [0, 1, 2d, 3d^2]
        RowVector4d T;
        T << 0, 1, 2 * d, 3 * pow(d, 2);

        // 앞 스플라인의 1차 미분 - 뒤 스플라인의 1차 미분 = 0
        for (int k = 0; k < 4; ++k) {
            M(row, 4 * j + k) = T(k);
            M(row, 4 * (j + 1) + k) = -T(k);
        }
    }

    // 2차 미분 연속성 제약: 인접 스플라인 접점에서 2차 미분값이 같아야 함
    for (int j = 0; j < no_splines - 1; ++j) {
        int row = 3 * no_splines - 1 + j;
        double d = ds(j);

        // 2차 미분 다항식 벡터: p''(d) = [0, 0, 2, 6d]
        RowVector4d T;
        T << 0, 0, 2, 6 * d;

        // 앞 스플라인 2차 미분 - 뒤 스플라인 2차 미분 = 0
        for (int k = 0; k < 4; ++k) {
            M(row, 4 * j + k) = T(k);
            M(row, 4 * (j + 1) + k) = -T(k);
        }
    }

    // 시작점 방향 psi_s가 주어지면 방향 제약 추가 (1차 미분 방향 제약)
    if (!std::isnan(psi_s)) {
        M.row(4 * no_splines - 2).setZero();
        M(4 * no_splines - 2, 1) = 1;           // 첫 스플라인 1차항 계수만 사용
        b_x(4 * no_splines - 2) = cos(psi_s);  // 방향 벡터 x 성분
        b_y(4 * no_splines - 2) = sin(psi_s);  // 방향 벡터 y 성분
    }

    // 끝점 방향 psi_e가 주어지면 끝점 1차 미분 방향 제약 추가
    if (!std::isnan(psi_e)) {
        double d = ds(no_splines - 1);
        M.row(4 * no_splines - 1).setZero();

        // 마지막 스플라인 구간의 1차 미분 다항식 계수 적용
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 1) = 1;
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 2) = 2 * d;
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 3) = 3 * pow(d, 2);

        b_x(4 * no_splines - 1) = cos(psi_e);
        b_y(4 * no_splines - 1) = sin(psi_e);
    }

    // 연립방정식 M * coeffs = b_x, b_y 풀기 (Least squares 해법)
    VectorXd x_les = M.colPivHouseholderQr().solve(b_x);
    VectorXd y_les = M.colPivHouseholderQr().solve(b_y);

    // 1차원 벡터를 각 스플라인 구간별 4개 계수 행렬로 변환 (전치 필요)
    MatrixXd coeffs_x = Map<MatrixXd>(x_les.data(), 4, no_splines).transpose();
    MatrixXd coeffs_y = Map<MatrixXd>(y_les.data(), 4, no_splines).transpose();

    // 각 구간의 법선 벡터 계산: 법선 벡터 = (-dy, dx), dx, dy는 1차항 계수 (속도 성분)
    MatrixXd normvec(no_splines, 2);
    for (int i = 0; i < no_splines; ++i) {
        double dx = coeffs_x(i, 1);
        double dy = coeffs_y(i, 1);
        normvec(i, 0) = -dy;
        normvec(i, 1) = dx;
    }

    // 법선 벡터 크기 계산 및 정규화 (단위 벡터로 변환)
    VectorXd norms = normvec.rowwise().norm();
    MatrixXd normvec_normalized(no_splines, 2);
    for (int i = 0; i < no_splines; ++i) {
        normvec_normalized.row(i) = normvec.row(i) / norms(i);
    }

    // 결과 구조체 반환
    return {coeffs_x, coeffs_y, M, normvec_normalized};
}

int main() {
    // 예시 경로 4점 (x,y)
    MatrixXd path(4, 2);
    path << 0, 0,
            1, 2,
            3, 3,
            4, 0;

    // 스플라인 계수 계산 함수 호출
    SplineResult result = calcSplines(path);

    // 결과 출력
    cout << "Coefficients X:\n" << result.coeffs_x << endl;
    cout << "Coefficients Y:\n" << result.coeffs_y << endl;
    cout << "Normalized Normals:\n" << result.normvec_normalized << endl;

    return 0;
}
