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

        // 2. 스케일링
        // ds라는 거리 벡터가 계산된 상태에서, 각 구간 간의 길이 비율(scaling)을 계산 해 VectorXd scaling에 저장
        // scaling 벡터 크기 : no_splines - 1 (두 스플라인 사이의 연결점 개수 == 스플라인 개수 - 1)

        // 이후 스플라인 연립방정식 구성 시 사용 (연속성 조건에 반영)
        VectorXd scaling = VectorXd::Ones(no_splines - 1);  // 기본값: 모두 1.0
        if (use_dist_scaling) {
            for (int i = 0; i < no_splines - 1; ++i) {
                scaling(i) = ds(i) / ds(i + 1);  // 인접 구간 간의 거리 비율 계산
            }
        }

        // 3. 행렬 M 구성
        // 각 스플라인은 x, y 각각에 대해 4개의 계수(a0 ~ a3) 가짐
        // 따라서 N개의 점이면, N-1개의 스플라인 구간 * 4개의 계수 = 총 4(N-1)개의 계수 발생

        // 위치 조건(각 스플라인은 시작점과 끝점에서 위치가 정확히 맞아야 함) : 식 2*(N-1)개
        // 1차 도함수 연속(기울기 연속) : 식 (N-2)개
        // 2차 도함수 연속(곡률 연속) : 식 (N-2)개
        // 양 끝 점(경계 조건) : 식 2개
    

        SplineResult result;
        result.coeffs_x = MatrixXd::Zero(no_splines, 4);
        result.coeffs_y = MatrixXd::Zero(no_splines, 4);
        result.M = MatrixXd::Zero(1, 1); // dummy
        result.normvec_normalized = MatrixXd::Zero(no_splines, 2);

    return result;

}


int main() {
    using namespace Eigen;

    // 간단한 테스트 경로: 삼각형 형태
    MatrixXd path(3, 2);
    path << 0, 0,
            1, 1,
            0, 2;

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
