#include "graph_planner.hpp"
#include "config.h"

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
static void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale = 0.5) {
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
void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params) {
    plt::clf();

    // 트랙 경계
    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // 레이싱 라인 + 샘플링된 포인트 + 헤딩
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}, {"label", "Raceline"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}, {"label", "Sampled Raceline"}});
    plotHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);

    // 노드
    plotHeading(nodesPerLayer);

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