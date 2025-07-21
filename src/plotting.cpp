#include "graph_planner.hpp"
#include <unordered_set>

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
    int layer_idx, node_idx;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            double dx = scale * cos(node.psi + M_PI_2);
            double dy = scale * sin(node.psi + M_PI_2);

            node_x.push_back(node.x);
            node_y.push_back(node.y);
            plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}});

            ostringstream label;
            label.precision(2);
            label << fixed << layer_idx << ", " << node_idx;

            plt::text(layer_idx, node_idx, label.str());
            

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
            layer_idx++;
    }
}

vector<Vector2d> generateSplinePoints(const SplineResult& spline, int num_points = 20) {
    vector<Vector2d> points;
    points.reserve(num_points + 1);

    for (int i = 0; i <= num_points; ++i) {
        double t = static_cast<double>(i) / num_points;  

        double t2 = t * t;
        double t3 = t2 * t;

        double x = spline.coeffs_x(0, 0) + spline.coeffs_x(0, 1) * t + spline.coeffs_x(0, 2) * t2 + spline.coeffs_x(0, 3) * t3;
        double y = spline.coeffs_y(0, 0) + spline.coeffs_y(0, 1) * t + spline.coeffs_y(0, 2) * t2 + spline.coeffs_y(0, 3) * t3;

        points.emplace_back(x, y);
    }
    
    return points;
}

void plotAllSplines(const SplineMap& splineMap, const NodeMap& nodesPerLayer, const string &color) {
    unordered_set<int> plotted_layers;
    for (const auto& [edge_key, spline] : splineMap) {
        const auto& [startKey, endKey] = edge_key;  
        int start_layer = startKey.first;
        int start_idx = startKey.second;
        int end_layer = endKey.first;
        int end_idx = endKey.second;

        const Node& startNode = nodesPerLayer.at(start_layer).at(start_idx);
        const Node& endNode = nodesPerLayer.at(end_layer).at(end_idx);

        vector<Vector2d> spline_points = generateSplinePoints(spline);

        DVector xs, ys;
        for (const auto& pt : spline_points) {
            xs.push_back(pt.x());
            ys.push_back(pt.y());
        }

        // 선 그리기
        plt::plot(xs, ys, {{"color", color}});

        // 레이어당 한 번만 텍스트 표시
        if (plotted_layers.find(start_layer) == plotted_layers.end() && !spline_points.empty()) {
            const auto& mid_pt = spline_points[spline_points.size() / 2];
            string layer_label = "L" + std::to_string(start_layer);
            plt::text(mid_pt.x(), mid_pt.y(), layer_label);
            plotted_layers.insert(start_layer);
        }
    }

    plt::title("All Spline Paths");
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::axis("equal");
    plt::show();
}


void visual(const Graph& edgeList, const NodeMap& nodesPerLayer, const SplineMap& splineMap, const string &color) {
    plt::clf();

    DVector x_bound_l = sampling_map[__x_bound_l];
    DVector y_bound_l = sampling_map[__y_bound_l];

    DVector x_bound_r = sampling_map[__x_bound_r];
    DVector y_bound_r = sampling_map[__y_bound_r];
    // 첫 점을 맨 뒤에 추가
    x_bound_l.push_back(x_bound_l.front());
    y_bound_l.push_back(y_bound_l.front());
    
    x_bound_r.push_back(x_bound_r.front());
    y_bound_r.push_back(y_bound_r.front());

    plt::plot(x_bound_l, y_bound_l, {{"color", "orange"}});
    plt::plot(x_bound_r, y_bound_r, {{"color", "orange"}});

    plt::plot(sampling_map[__x_bound_l], sampling_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(sampling_map[__x_bound_r], sampling_map[__y_bound_r], {{"color", "orange"}});
    
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
    plotHeading(nodesPerLayer);

    plotAllSplines(splineMap, nodesPerLayer, color);

    // plt::title("Track");
    // plt::grid(true);
    // plt::axis("equal");
    // plt::show();  
}