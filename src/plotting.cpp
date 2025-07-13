#include "graph_planner.hpp"

void plotHeading(const DVector &x,
                 const DVector &y,
                 const DVector &psi,
                 double scale = 0.5)
{
    double dx, dy;
    double theta, arrow_len;
    double angle;
    double x_arrow1, y_arrow1;
    double x_arrow2, y_arrow2;

    for (size_t i = 0; i < x.size(); ++i) {
        dx = scale * cos(psi[i] + M_PI_2);
        dy = scale * sin(psi[i] + M_PI_2);

        // psi 방향 
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        // plt::plot(x_line, y_line, {{"color", "green"}});

        #if 0
        // 화살촉 
        theta = atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;  // 30 degrees

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);

        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        // 화살촉 그리기 
        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "green"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "green"}});
        #endif

    }
        // raceline 좌표와 psi 프린팅 
        #if 1
        for (size_t i = 0; i < x.size(); ++i) {
        ostringstream label;
        label.precision(2);
        label << fixed << "(" << x[i] << ", " << y[i] << ")\nψ=" << psi[i];

        plt::text(x[i], y[i], label.str());
    }
        #endif
}

void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector x_line, y_line;
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            double dx = scale * cos(node.psi + M_PI_2);
            double dy = scale * sin(node.psi + M_PI_2);

            node_x.push_back(node.x);
            node_y.push_back(node.y);
            plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}});

            #if 0
            x_line = {node.x, node.x + dx};
            y_line = {node.y, node.y + dy};
            plt::plot(x_line, y_line, {{"color", "purple"}});

            // 화살촉 (arrowhead)
            
            double theta = atan2(dy, dx);
            double arrow_len = 0.2 * scale;
            double angle = M_PI / 6.0;

            double x_arrow1 = node.x + dx - arrow_len * cos(theta - angle);
            double y_arrow1 = node.y + dy - arrow_len * sin(theta - angle);

            double x_arrow2 = node.x + dx - arrow_len * cos(theta + angle);
            double y_arrow2 = node.y + dy - arrow_len * sin(theta + angle);

            plt::plot({node.x + dx, x_arrow1}, {node.y + dy, y_arrow1}, {{"color", "purple"}});
            plt::plot({node.x + dx, x_arrow2}, {node.y + dy, y_arrow2}, {{"color", "purple"}});
            #endif
        }
        
    }
}

void plotSplinesFromMap(const SplineMap& splineMap, const NodeMap& nodesPerLayer, int num_samples = 20) {
    std::vector<double> spline_x_pts;
    std::vector<double> spline_y_pts;

    for (const auto& [edgeKey, splineRes] : splineMap) {
        spline_x_pts.clear();
        spline_y_pts.clear();

        const IPair& startKey = edgeKey.first;
        const IPair& endKey = edgeKey.second;

        const int start_layer = startKey.first;
        const int start_idx   = startKey.second;
        const int end_layer   = endKey.first;
        const int end_idx     = endKey.second;

        // 유효성 검사
        if (start_layer >= nodesPerLayer.size() || end_layer >= nodesPerLayer.size()) continue;
        if (start_idx >= nodesPerLayer[start_layer].size() || end_idx >= nodesPerLayer[end_layer].size()) continue;

        const Node& startNode = nodesPerLayer[start_layer][start_idx];
        const Node& endNode   = nodesPerLayer[end_layer][end_idx];

        const auto& coeffs_x = splineRes.coeffs_x;
        const auto& coeffs_y = splineRes.coeffs_y;

        // 일반적으로 1개 구간만 있다고 가정 (row 0)
        for (int i = 0; i <= num_samples; ++i) {
            double t = static_cast<double>(i) / num_samples;

            double x = coeffs_x(0, 0)
                     + coeffs_x(0, 1) * t
                     + coeffs_x(0, 2) * t * t
                     + coeffs_x(0, 3) * t * t * t;

            double y = coeffs_y(0, 0)
                     + coeffs_y(0, 1) * t
                     + coeffs_y(0, 2) * t * t
                     + coeffs_y(0, 3) * t * t * t;

            spline_x_pts.push_back(x);
            spline_y_pts.push_back(y);
        }

        // 그리기
        plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}});
    }
}


void visual(const Graph& edgeList, const NodeMap& nodesPerLayer, const SplineMap& splineMap) {
    plt::clf();

    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // plt::plot(gtpl_map[__x_ref], gtpl_map[__y_ref], {{"color", "blue"}});
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}});

    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}});

    // plotHeading(sampling_map[__x_raceline],
    //             sampling_map[__y_raceline],
    //             sampling_map[__psi]);

    // plotHeading(sampling_map[__x_bound_l],
    //             sampling_map[__y_bound_l],
    //             psi_bound_l);

    // plotHeading(sampling_map[__x_bound_r],
    //             sampling_map[__y_bound_r],
    //             psi_bound_r);


    // 노드마다 psi확인할 수 있는 용도 
    // plotHeading(nodesPerLayer);
    plotSplinesFromMap(splineMap, nodesPerLayer);   

    plt::title("Track");
    plt::grid(true);
    plt::axis("equal");
    plt::show();  
}