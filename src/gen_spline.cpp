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
        

        // [DEBUG] M, b_x, b_y 출력
        cout << "[DEBUG] M 행렬:\n" << M << "\n";
        cout << "[DEBUG] b_x 벡터:\n" << b_x.transpose() << "\n";
        cout << "[DEBUG] b_y 벡터:\n" << b_y.transpose() << "\n";

        // 연립방정식 풀기
        VectorXd coeffs_x = M.colPivHouseholderQr().solve(b_x);
        VectorXd coeffs_y = M.colPivHouseholderQr().solve(b_y);

        // 결과 구조체 초기화
        SplineResult result;
        result.coeffs_x = MatrixXd::Zero(no_splines, 4);
        result.coeffs_y = MatrixXd::Zero(no_splines, 4);
        result.M = M;
        result.normvec_normalized = MatrixXd::Zero(no_splines, 2); // 아직 법선 벡터 미구현

        // 계수 복사
        for (int i = 0; i < no_splines; i++) {
            for (int j = 0; j < 4; j++) {
                result.coeffs_x(i, j) = coeffs_x(i * 4 + j);
                result.coeffs_y(i, j) = coeffs_y(i * 4 + j);
            }
        }

        return result;

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
