#include "graph_planner.hpp"

DMap gtpl_map;
DMap sampling_map;

void samplePointsFromRaceline(const DVector& kappa,
                              const DVector& dist,
                              double d_curve,
                              double d_straight,
                              double curve_th,
                              IVector& idx_array) {
    // idx_array: sampling된 raceline 위 인덱스

    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    for (size_t i = 0; i < n; ++i) {

        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist;
        }
        // cout << fabs(kappa[i]) << endl;
        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) {
            idx_array.push_back(static_cast<int>(i));

            if (fabs(kappa[i]) < curve_th) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist[i];
    }

    // for (size_t i=0; i < idx_array.size(); ++i) 
    //     cout << idx_array[i] << endl;
    // cout << "size: " << idx_array.size() << endl;
}

void genNode(NodeMap& nodesPerLayer,
            IVector& raceline_index_array,
            const double veh_width,
            float lat_resolution) {
    
    const int N = sampling_map[__alpha].size();
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);
        // cout << i << "번째 layer 길이" << (sampling_map[__width_left][i] + sampling_map[__width_right][i])<< endl;
        // cout << "layer 내에서 raceline index:" << raceline_index << endl;
        // cout << "-----" << endl;

        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);    // 기준선에서의 위치
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);  // 기준선에서 수직한 노멀 벡터 따라 노드 배치
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지
        int node_idx = 0;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution;  // num_nodes : 좌우 총 가능한 노드 수
        
        nodesPerLayer[i].resize(num_nodes); 
        
        // cout << i << "번째 layer의 node 개수는 " << num_nodes << endl;  
        // node별 loop 
        for (int idx = 0; idx < num_nodes; ++idx) {
            double alpha = start_alpha + idx * lat_resolution;
            // cout << idx << "번째 노드" << endl;
            // node의 좌표 계산.
            node_pos = ref_xy + alpha * norm_vec;

            node.x = node_pos.x();
            node.y = node_pos.y();      
            node.raceline = (node_idx == raceline_index);
            
            #if 0 
            if (idx == num_nodes - 1) {
                cout << i << "번째 레이어의 " << idx << "번째 노드" << endl;
                cout << sampling_map[__x_bound_r][idx] - alpha << endl;
            }
            #endif

            // psi 재계산  
            double psi_interp;
            if (node_idx < raceline_index) {
                
                if (abs(sampling_map[__psi_bound_l][i] - sampling_map[__psi][i]) >= M_PI) 
                {   
                    double bl = sampling_map[__psi_bound_l][i] + 2 * M_PI * (sampling_map[__psi_bound_l][i] < 0);
                    double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
                    psi_interp = bl + (p - bl) * node_idx / raceline_index;          
                }
                else {
                    psi_interp = sampling_map[__psi_bound_l][i] + (sampling_map[__psi][i] - sampling_map[__psi_bound_l][i]) * (node_idx+1) / raceline_index;
                }
                node.psi = normalizeAngle(psi_interp);
            }
            else if (node_idx == raceline_index) {
                psi_interp = sampling_map[__psi][i];
                node.psi = psi_interp;
            }
            else {
                int remain = num_nodes - raceline_index - 1;
                double t = static_cast<double>(node_idx - raceline_index) / max(remain, 1);  // 0 ~ 1
                psi_interp = sampling_map[__psi][i] + t * (sampling_map[__psi_bound_r][i] - sampling_map[__psi][i]);
                node.psi = normalizeAngle(psi_interp);
            }
            // cout << i << "번째 레이어의" <<node_idx << "번째 노드의 psi는" << node.psi << endl;

            nodesPerLayer[i][node_idx] = node;
            ++node_idx;

        }

    }
        // cout << i << "번째 Layer의" << endl;
        // for (size_t i =0; i < node_pos.size(); ++i) {        
        //     cout << i << "번째 Node" << endl;
        //     cout << node_pos[i] << endl;
        // }

}

void calcOfflineCost(SplineMap& splineMap,
                   IVector& raceline_index_array,
                   float w_curv_avg,
                   float w_curv_peak, 
                   float w_length, 
                   float lat_resolution,
                   float w_raceline, 
                   float w_raceline_sat) {
    if (splineMap.size() <= 0) {
        throw invalid_argument("SplineMap's Size is zero!!");
    }

    for (auto& [startPoint, endPoints] : splineMap) {
        for (auto& [endPoint, spline] : endPoints) {
            double offline_cost = 0.0;
            int end_layer = endPoint.first;
            int end_node = endPoint.second;

            // 디버깅용 
            // cout << "kappa: ";
            // for (int i = 0; i < spline.kappa.size(); ++i) cout << spline.kappa[i] << " ";
            // cout << endl;

            if (end_layer < 0 || end_layer >= raceline_index_array.size())
            {
                continue;
            }

            if (spline.kappa.size() == 0)
            {
                // cerr << "Empty kappa in spline!" << endl;
                continue;
            }

            double abs_kappa = spline.kappa.array().abs().sum();
            double s_length = spline.el_lengths.sum();
            // cout << "s_length: " << s_length << endl;
            // average curvature
            offline_cost += w_curv_avg * pow(abs_kappa / float(spline.kappa.size()), 2) * s_length;
            // peak curvature
            double max_min = std::abs(spline.kappa.array().maxCoeff() - spline.kappa.array().minCoeff());
            offline_cost += w_curv_peak * pow(max_min, 2) * s_length;

            // path length
            offline_cost += w_length * s_length;

            // raceline cost

            double raceline_dist = std::abs(raceline_index_array[end_layer] - end_node) * lat_resolution;
            double raceline_cost = min(w_raceline * s_length * raceline_dist, w_raceline_sat * s_length);

            offline_cost += raceline_cost;

            spline.cost = offline_cost;
            // cout << "(" << startPoint.first << ", " << startPoint.second << ") " << " -> " << "(" << end_layer << ", " << end_node << "): " << offline_cost << endl;
        }   
    } 
}

void getClosestNodes(const NodeMap& nodesPerLayer, IPair& closest_idx, const Vector2d& pos, int limit=1) {

    int num_nodes = 0;
    for (const auto& layer : nodesPerLayer) {
        num_nodes += layer.size();
    }

    MatrixXd node_xy(num_nodes, 2);
    int idx = 0;

    for (size_t i = 0; i < nodesPerLayer.size(); ++i) {
        for (size_t j = 0; j < nodesPerLayer[i].size(); ++j) {
            const Node& node = nodesPerLayer[i][j];
            node_xy(idx, 0) = node.x;
            node_xy(idx, 1) = node.y;
            ++idx;
        }
    }   
    // pos(2, 1) -> pos.transpose() -> (1, 2)
    MatrixXd diff = node_xy.rowwise() - pos.transpose();
    VectorXd dist2 = diff.rowwise().squaredNorm();
    vector<tuple<double, int, int>> dist_info;

    int re_idx = 0;
    for (size_t i = 0; i < nodesPerLayer.size(); ++i) {
        for (size_t j = 0; j < nodesPerLayer[i].size(); ++j) {
            dist_info.emplace_back(dist2(re_idx++), i, j);
        }
    }

    // 최소 거리 limit개만 앞으로 정렬
    nth_element(dist_info.begin(), dist_info.begin() + limit, dist_info.end());

    // 결과 저장
    for (int k = 0; k < limit; ++k) {
        auto [dist, i, j] = dist_info[k];
        closest_idx = make_pair(i, j);
        cout << "Closest node: layer=" << i << ", idx=" << j << endl;
    }
}

void setInitialPos(const NodeMap &nodesPerLayer,
                   const IVector &raceline_index_array, 
                   const float &max_heading_offset, 
                   const float &stepsize_approx,
                   float &veh_width)
{
    // 현재 pos, heading 
    double dx, dy;
        
    Vector2d initial_pos(sampling_map[__x_raceline][1], sampling_map[__y_raceline][1]);
    float vel_est = 0.0;

    // set start pos 
    if (!checkInsideBounds(initial_pos, veh_width)) {
        throw out_of_range("start pos is not in bounds");
    }
    
    IPair closest_idx;
    getClosestNodes(nodesPerLayer, closest_idx, initial_pos);
    int start_layer = closest_idx.first;
    int start_node = closest_idx.second;
    double start_heading = nodesPerLayer[closest_idx.first][closest_idx.second].psi;

    int end_layer = (closest_idx.first + 2) % (nodesPerLayer.size() - 1);

    for (int layer_idx = start_layer; layer_idx < end_layer;++layer_idx) {
        if (layer_idx != start_layer) {
            start_node = raceline_index_array[layer_idx];
            start_heading = nodesPerLayer[layer_idx][start_node].psi;
        }

        int goal_layer = layer_idx + 1;
        int goal_node = raceline_index_array[goal_layer];
        double goal_heading = nodesPerLayer[goal_layer][goal_node].psi;
        double heading_diff = std::abs(start_heading - goal_heading);
        Vector2d start_pos(nodesPerLayer[layer_idx][start_node].x, nodesPerLayer[layer_idx][start_node].y);
        Vector2d end_pos(nodesPerLayer[goal_layer][goal_node].x, nodesPerLayer[goal_layer][goal_node].y);

        if (heading_diff > M_PI) {
            heading_diff = std::abs(2*M_PI - heading_diff);
            // cout <<"hi" << endl;
        }
        if (heading_diff > max_heading_offset) {
            // cout << "heading_diff: " << heading_diff << ", " << max_heading_offset << endl;
            // cerr << "Heading mismatch between vehicle and track grid, check if vehicle oriented correctly!" << endl;
        }

        MatrixXd path(2, 2);
        path.block<1, 2>(0, 0) = start_pos.transpose();
        path.block<1, 2>(1, 0) = end_pos.transpose();

        auto result = calcSplines(path, start_heading, goal_heading); // coeffs_x, coeffs_y, kappa, el_lengths, cost
        auto [kappa, psi] = interpSplines(result->coeffs_x, result->coeffs_y, stepsize_approx, veh_width);

        // kappa와 psi의 크기 확인
        if (kappa.size() == 0) {
            cerr << "Error: kappa is empty!" << endl;
        }

        if (psi.size() == 0) {
            cerr << "Error: psi is empty!" << endl;
            return;
        }
        plotSpline(*result, "blue");    
    }
    cout << "Goal node: layer=" << end_layer << ", idx= " << raceline_index_array[end_layer] << endl;
    #if 0
    ActionSet actionSet;
    actionSet.action_id = "straight";

    MatrixXd coeffs_all(2, 4);
    coeffs_all.row(0) = result->coeffs_x;
    coeffs_all.row(1) = result->coeffs_y;

    // cout << coeffs_all.row(0) << endl;

    actionSet.coeffs.push_back(coeffs_all);
    // cout << result->el_lengths.rows() << result->el_lengths.cols();
    VectorXd el_lengths_all(2);
    el_lengths_all(0) = result->el_lengths(0); // 거리
    el_lengths_all(1) = 0.0;

    MatrixXd action_param(1, 5);
    action_param(0, 0) = start_pos.x();  
    action_param(0, 1) = start_pos.y();   
    action_param(0, 2) = psi(0);             
    action_param(0, 3) = kappa(0);           
    action_param(0, 4) = el_lengths_all(0);
    #endif

}

int main() {
    clock_t s_time, f_time;
    s_time = clock();

    IVector idx_sampling;
    Offline_Params params;

    unique_ptr<string> track = Load("include/driving_task.ini");

    // 3. 자동 경로 설정
    string map_file_in  = "inputs/traj_ltpl_cl_" + *track + ".csv";
    string map_file_out = "outputs/" + *track + "_out.csv";

    // global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    readDMapFromCSV(map_file_in, gtpl_map);

    addDVectorToMap(gtpl_map, "bound_r");
    addDVectorToMap(gtpl_map, "bound_l");
    addDVectorToMap(gtpl_map, "raceline");
    addDVectorToMap(gtpl_map, "delta_s");

    writeDMapToCSV(map_file_out, gtpl_map);
    
    // layer 간격을 위한 raceline points sampling 
    samplePointsFromRaceline(gtpl_map[__kappa],
                             gtpl_map[__delta_s],
                             params.LON_CURVE_STEP,
                             params.LON_STRAIGHT_STEP,
                             params.CURVE_THR,
                             idx_sampling);

    // cout << "idx size:" << idx_sampling.size() << endl;
    
    for (const auto& [key, vec] : gtpl_map) {
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < vec.size()) {
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }
    // writeDMapToCSV("inputs/sampling_map", sampling_map);
    // map_size(sampling_map); // (51, 3)

    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);
    
    // map_size(sampling_map); // (51, 4)

    // 추후 저장될 예정 
    calcHeading(sampling_map[__x_raceline],
                sampling_map[__y_raceline],
                sampling_map[__psi]);
    
    // 여기서 계산되는 sampling된 bound_l, r은 node 생성 시에만 쓰인다. 
    calcHeading(sampling_map[__x_bound_l],
                sampling_map[__y_bound_l],
                sampling_map[__psi_bound_l]);

    calcHeading(sampling_map[__x_bound_r],
                sampling_map[__y_bound_r],
                sampling_map[__psi_bound_r]);  

    // sampling_map[__psi_bound_l] = psi_bound_l;
    // sampling_map[__psi_bound_r] = psi_bound_r;

    NodeMap nodesPerLayer;
    IVector raceline_index_array;
    // Vector2d node_pos;

    genNode(nodesPerLayer,
            raceline_index_array,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    // for (auto index : raceline_index_array) {
    //     cout << index << endl;
    // }

    // sampling points' info 
    // writeDMapToCSV("inputs/sampling_map.csv", sampling_map);

    Graph graph_wp; // a graph of waypoints
    SplineMap splineMap;
    genEdges(nodesPerLayer,
             graph_wp,
             splineMap,
             raceline_index_array,
             params.VEH_WIDTH,
             params.LAT_OFFSET,
             params.LAT_RESOLUTION,
             params.CURVE_THR,
             params.MAX_LAT_STEPS,
             params.STEPSIZE_APPROX,
             params.MIN_VEL_RACE,
             params.MAX_LATERAL_ACCEL,
             params.VEH_TURN);

    calcOfflineCost(splineMap,
                   raceline_index_array,
                   params.W_CURV_AVG,
                   params.W_CURV_PEAK, 
                   params.W_LENGTH,
                   params.LAT_RESOLUTION, 
                   params.W_RACELINE, 
                   params.W_RACELINE_SAT);

    setInitialPos(nodesPerLayer,
                  raceline_index_array, 
                  params.MAX_HEADING_OFFSET, 
                  params.STEPSIZE_APPROX, 
                  params.VEH_WIDTH);

    f_time = clock();

    // graph_wp.printGraph();

    // visual process 
    cout << (double)(f_time - s_time) / CLOCKS_PER_SEC << "s 소요" << endl;
    
    // printSplineInfo(splineMap, nodesPerLayer);

    visual(graph_wp, nodesPerLayer, splineMap, "gray");

    return 0;
}