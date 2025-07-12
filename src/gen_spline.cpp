#include "graph_planner.hpp"
#include "config.h"

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

        // 1. 연립방정식 풀기 (QR 분해로 안정적)
        VectorXd x_flat = M.colPivHouseholderQr().solve(b_x);  // (4 * N x 1)
        VectorXd y_flat = M.colPivHouseholderQr().solve(b_y);  // (4 * N x 1)

        // 2. 계수 벡터를 (N x 4) 행렬로 reshape
        MatrixXd coeffs_x = Map<Matrix<double, 4, Dynamic>>(x_flat.data(), 4, no_splines).transpose();
        MatrixXd coeffs_y = Map<Matrix<double, 4, Dynamic>>(y_flat.data(), 4, no_splines).transpose();

        // 3. 법선 벡터 계산 (노멀 벡터)
        MatrixXd normvec(no_splines, 2);
        for (int i = 0; i < no_splines; ++i) {
            double dx = coeffs_x(i, 1);  // a₁: x방향 도함수 시작점
            double dy = coeffs_y(i, 1);  // a₁: y방향 도함수 시작점
            normvec(i, 0) = -dy;
            normvec(i, 1) = dx;
        }

        VectorXd norms = normvec.rowwise().norm();
        MatrixXd normvec_normalized = normvec.array().colwise() / norms.array();

        // 4. 결과 구조체에 저장
        SplineResult result;
        result.coeffs_x = coeffs_x;
        result.coeffs_y = coeffs_y;
        result.M = M;
        result.normvec_normalized = normvec_normalized;


        return result;

}

void genEdges() {
    
}


int main() {
    using namespace Eigen;

    // 간단한 테스트 경로: 삼각형 형태
    MatrixXd path(5, 2);
    path << 0, 0,
            1, 2,
            3, 3,
            6, 2,
            10, 0;
    

    // 시작, 끝 방향 (라디안) – 테스트용 값
    double psi_s = M_PI / 2.0;
    double psi_e = M_PI / 2.0;

    // calcSplines 호출 (el_lengths 없이 테스트)
    SplineResult result = calcSplines(path, nullptr, psi_s, psi_e, true);

    // 결과 출력
    std::cout << "X 계수:\n" << result.coeffs_x << "\n\n";
    std::cout << "Y 계수:\n" << result.coeffs_y << "\n\n";
    std::cout << "M 행렬 크기: " << result.M.rows() << "x" << result.M.cols() << "\n";
    std::cout << "정규화된 노멀 벡터:\n" << result.normvec_normalized << "\n";

    return 0;
}
