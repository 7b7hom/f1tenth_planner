#include "graph_planner.hpp"
#include <unordered_set>

void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale = 0.5) {
    double dx, dy;
    double theta, arrow_len;
    double angle;
    double x_arrow1, y_arrow1;
    double x_arrow2, y_arrow2;

    for (size_t i = 0; i < x.size(); ++i) {
        dx = scale * cos(psi[i] + M_PI_2);
        dy = scale * sin(psi[i] + M_PI_2);
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        plt::plot(x_line, y_line, {{"color", "green"}});

        #if 1 // 화살촉 그리기
        theta = atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);
        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "green"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "green"}});
        #endif
    }
}

// NodeMap에 저장된 모든 노드들을 보라색 점으로 플로팅하고, 각 노드의 헤딩을 화살표로 시각화
void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            node_x.push_back(node.x);
            node_y.push_back(node.y);
        }
    }
    plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}, {"label", "Nodes"}});
}

void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params) {
    plt::clf();

    // 트랙 경계선
    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // 레이싱 라인 및 샘플링된 포인트
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}, {"label", "Raceline"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}, {"label", "Sampled Raceline"}});
    plotHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);

    plotHeading(nodesPerLayer);
    
    DVector spline_x_pts; 
    DVector spline_y_pts;

    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& current_node : layer_nodes) {
            IPair src_key = std::make_pair(current_node.layer_idx, current_node.node_idx);
            IPairVector child_nodes_idx; 

            try {
                graph.getChildNodes(src_key, child_nodes_idx);
            } catch (const std::runtime_error& e) {
                continue;
            }
            
            for (const auto& dst_key : child_nodes_idx) {
                spline_x_pts.clear(); 
                spline_y_pts.clear(); 

                size_t next_layer_idx = (current_node.layer_idx + 1) % nodesPerLayer.size();

                int dest_node_idx = dst_key.second;
                
                if (dest_node_idx < 0 || dest_node_idx >= static_cast<int>(nodesPerLayer[next_layer_idx].size())) {
                    continue;
                }

                const Node& next_node = nodesPerLayer[next_layer_idx][dest_node_idx];

                // 스플라인 경로 생성
                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y,
                            next_node.x, next_node.y;
                
                SplineResult res;
                try {
                    res = calcSplines(current_node, next_node);
                } catch (const std::exception& e) {
                    continue;
                }

                const int num_spline_segments = 10; 
                for (int k = 0; k <= num_spline_segments; ++k) {
                    double t_eval = static_cast<double>(k) / num_spline_segments;
                    Vector2d pos = computeSplinePosition(res.coeffs_x.row(0), res.coeffs_y.row(0), t_eval);
                    spline_x_pts.push_back(pos.x());
                    spline_y_pts.push_back(pos.y());
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