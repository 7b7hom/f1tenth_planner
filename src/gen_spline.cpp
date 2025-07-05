#include "graph_planner.hpp"
#include "config.h"

// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;          // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;          // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                 // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;// 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

// path : x,y 좌표 배열
// el_lengths 전달 안 하면 nullptr (점 사이 거리는 필수 x)
// psi_s, psi_e도 전달 안 하면 NAN (필수 x)
SplineResult calcSplines(const MatrixXd& path, const VectorXd* el_lengths = nullptr,
                         double psi_s = NAN, double psi_e = NAN,
                         bool use_dist_scaling = true) {

        // 폐곡선 여부 판단
        MatrixXd updated_path = path;

        bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6 && std::isnan(psi_s);

        if (closed) {
            // 폐곡선 유지: 시작점 다시 한 번 넣기
            updated_path.conservativeResize(path.rows() + 1, Eigen::NoChange);
            updated_path.row(updated_path.rows() - 1) = path.row(0);
        }

        int no_splines = updated_path.rows() - 1;



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
