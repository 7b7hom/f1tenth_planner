#include "graph_planner.hpp"
#include "config_modena.h"

extern DMap gtpl_map;
extern DMap sampling_map;

struct SplineResult;
struct SplinePoint;
SplineResult calcSplines(const MatrixXd& path,
                         const VectorXd* el_lengths_ptr,
                         double psi_s, double psi_e, bool use_dist_scaling);
SplinePoint evaluateSpline(const RowVector4d& coeff_x,
                           const RowVector4d& coeff_y,
                           double t, double ds_current, bool normalized_t);

// ---- 내부 유틸: 헤딩 화살표 ----
static void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale = 0.1) {
    for (size_t i = 0; i < x.size(); ++i) {
        double dx = scale * cos(psi[i] + M_PI_2);
        double dy = scale * sin(psi[i] + M_PI_2);
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        plt::plot(x_line, y_line, {{"color", "red"}});

        double theta = atan2(dy, dx);
        double arrow_len = 0.2 * scale;
        double angle = M_PI / 6.0;

        double x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        double y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);
        double x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        double y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "red"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "red"}});
    }
}

// ---- 내부 유틸: 노드 점 찍기 ----
static void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            node_x.push_back(node.x);
            node_y.push_back(node.y);
        }
    }
    plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}, {"label", "Nodes"}});
}

// ---- 공개 함수: 전체 시각화 ----
void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params,
            const Vector2d& pos_est, double heading_est, const ITuple& start_key, int next_idx) {
    plt::clf();

    // 트랙 경계
    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // 레이싱 라인 + 샘플링된 포인트 + 헤딩
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}, {"label", "Raceline"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}, {"label", "Sampled Raceline"}});
    plotHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);

    // 노드
    //plotHeading(nodesPerLayer);

    // --- 노드 점과 노드 헤딩을 그리는 루프 추가 ---
    DVector node_x, node_y, node_psi;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            node_x.push_back(node.x);
            node_y.push_back(node.y);
            node_psi.push_back(node.psi);
        }
    }
    // 노드 점 그리기
    plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}, {"label", "Nodes"}});
    // 노드 헤딩 그리기 (scale을 0.5로 설정하여 레이스라인보다 짧게)
    plotHeading(node_x, node_y, node_psi, 1.0);

    // --- startpos 시각화 (하이라이트) ---
    // start_key의 layer 인덱스가 0 이상일 때만 (즉, set_startpos가 성공했을 때만) 실행
    if (get<0>(start_key) >= 0) {
        // startpos 관련 모든 시각화 요소를 이 블록 안에 넣습니다.
        
        // 1. startpos 추정 위치 시각화 (파란색)
        plt::scatter(std::vector<double>{pos_est.x()}, std::vector<double>{pos_est.y()}, 100.0, {{"color", "blue"}, {"label", "Start Pos Est"}});
        
        // 2. startpos 추정 헤딩 시각화
        double heading_arrow_len = 2.0;
        double dx = heading_arrow_len * cos(heading_est + M_PI_2);
        double dy = heading_arrow_len * sin(heading_est + M_PI_2);
        plt::plot(std::vector<double>{pos_est.x(), pos_est.x() + dx}, std::vector<double>{pos_est.y(), pos_est.y() + dy}, {{"color", "blue"}, {"linewidth", "2"}});
        
        int start_layer_idx = get<0>(start_key);
        int start_node_idx = get<1>(start_key);

        // 3. start_key 노드 하이라이트 (라임색)
        if (start_layer_idx < nodesPerLayer.size() && start_node_idx < nodesPerLayer[start_layer_idx].size()) {
            const Node& start_node = nodesPerLayer[start_layer_idx][start_node_idx];
            plt::scatter(std::vector<double>{start_node.x}, std::vector<double>{start_node.y}, 150.0, {{"color", "lime"}, {"marker", "o"}, {"label", "Start Node"}});

            // 4. next_idx 노드 하이라이트 (노란색)
            int next_layer_idx = (start_layer_idx + 1) % nodesPerLayer.size();
            if (next_layer_idx < nodesPerLayer.size() && next_idx >= 0 && next_idx < nodesPerLayer[next_layer_idx].size()) {
                const Node& next_node = nodesPerLayer[next_layer_idx][next_idx];
                plt::scatter(std::vector<double>{next_node.x}, std::vector<double>{next_node.y}, 150.0, {{"color", "yellow"}, {"marker", "o"}, {"label", "Next Node"}});

                // 5. 시작 노드에서 다음 노드까지의 경로 시각화 (자홍색)
                MatrixXd spline_path(2, 2);
                spline_path << start_node.x, start_node.y, next_node.x, next_node.y;
                VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();
                try {
                    SplineResult res = calcSplines(spline_path, &el_lengths, start_node.psi, next_node.psi, true);
                    DVector path_x, path_y;
                    const int num_seg = 10;
                    for (int k = 0; k <= num_seg; ++k) {
                        double t = double(k) / num_seg;
                        SplinePoint sp = evaluateSpline(res.coeffs_x.row(0), res.coeffs_y.row(0), t, res.ds(0), true);
                        path_x.push_back(sp.x);
                        path_y.push_back(sp.y);
                    }
                    plt::plot(path_x, path_y, {{"color", "magenta"}, {"linewidth", "3"}, {"label", "Start Path"}});
                } catch (...) {
                    // 스플라인 생성 실패 시 무시
                }
            }
        }
    }

    // 그래프 엣지(스플라인)
    DVector spline_x_pts, spline_y_pts;

    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& current_node : layer_nodes) {
            ITuple src_key(current_node.layer_idx, current_node.node_idx);
            IVector child_nodes_idx;
            try {
                graph.getChildIdx(src_key, child_nodes_idx);
            } catch (...) {
                continue;
            }

            for (int dest_node_idx : child_nodes_idx) {
                // 시작 경로와 겹치면 그리지 않음
                if (get<0>(start_key) == current_node.layer_idx && get<1>(start_key) == current_node.node_idx && next_idx == dest_node_idx) {
                    continue;
                }
                
                spline_x_pts.clear();
                spline_y_pts.clear();

                size_t next_layer_idx = (current_node.layer_idx + 1) % nodesPerLayer.size();
                if (dest_node_idx < 0 || dest_node_idx >= (int)nodesPerLayer[next_layer_idx].size()) continue;

                const Node& next_node = nodesPerLayer[next_layer_idx][dest_node_idx];

                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y,
                               next_node.x, next_node.y;

                double psi_s = current_node.psi;
                double psi_e = next_node.psi;

                Eigen::VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

                SplineResult res;
                try {
                    res = calcSplines(spline_path, &el_lengths, psi_s, psi_e, true);
                } catch (...) {
                    continue;
                }

                const int num_seg = 10;
                for (int k = 0; k <= num_seg; ++k) {
                    double t = double(k) / num_seg;
                    auto sp = evaluateSpline(res.coeffs_x.row(0), res.coeffs_y.row(0), t, res.ds(0), true);
                    spline_x_pts.push_back(sp.x);
                    spline_y_pts.push_back(sp.y);
                }
                plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}});
            }
        }
    }

    plt::title("Track and Planned Graph");
    plt::grid(true);
    plt::axis("equal");
    plt::legend();
    plt::show();
}