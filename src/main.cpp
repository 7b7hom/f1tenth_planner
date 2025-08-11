#include "graph.h"
#include "spline.h"
// #include "graph_planner.hpp"

void visual(DMap &gtMap, DMap &stMap, const NodeMap &nodesPerLayer, const SplineMap &splineMap);
void plotSpline(const Spline& spline, const string& color);

IVector samplePointsFromRaceline(const DVector& kappa,
                              const DVector& dist,
                              const YAML::Node& params) {

    IVector idx_sampling;
    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    // params
    float curve_thr = params["lattice"]["curve_thr"].as<float>();
    float d_straight = params["lattice"]["d_straight"].as<float>();
    float d_curve = params["lattice"]["d_curve"].as<float>();

    for (size_t i = 0; i < n; ++i) {
        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_thr) {
            next_dist = cur_dist;
        }
        // cout << fabs(kappa[i]) << endl;
        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) {
            idx_sampling.push_back(static_cast<int>(i));

            if (fabs(kappa[i]) < curve_thr) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist[i];
    }

    // for (size_t i=0; i < idx_sampling.size(); ++i) 
    //     cout << idx_sampling[i] << endl;
    // cout << "size: " << idx_sampling.size() << endl;

    return idx_sampling;
}


auto genNode(DMap &stMap, const YAML::Node &params) -> pair<NodeMap, IVector> {
    NodeMap nodesPerLayer;
    IVector raceline_index_array;
    
    const int N = stMap[NORM_L].size();
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    
    // params
    float veh_width = params["vehicle"]["veh_width"].as<float>();
    float lat_resolution = params["lattice"]["lat_resolution"].as<float>();
    
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((stMap[WIDTH_L][i] + stMap[NORM_L][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);

        Vector2d ref_xy(stMap[POS_X][i], stMap[POS_Y][i]);    // 기준선에서의 위치
        Vector2d norm_vec(stMap[NORM_X][i], stMap[NORM_Y][i]);  // 기준선에서 수직한 노멀 벡터 따라 노드 배치
        
        double start_alpha = stMap[NORM_L][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지
        int node_idx = 0;
        int num_nodes = (stMap[WIDTH_R][i] + stMap[WIDTH_L][i] - veh_width) / lat_resolution;  // num_nodes : 좌우 총 가능한 노드 수
        
        nodesPerLayer[i].resize(num_nodes); 
        
        // node별 loop 
        for (int idx = 0; idx < num_nodes; ++idx) {
            double alpha = start_alpha + idx * lat_resolution;
            // cout << idx << "번째 노드" << endl;
            // node의 좌표 계산.
            node_pos = ref_xy + alpha * norm_vec;

            node.x = node_pos.x();
            node.y = node_pos.y();      
            node.raceline = (node_idx == raceline_index);
            
            // psi 재계산  
            double psi_interp;
            if (node_idx < raceline_index) {
                
                if (abs(stMap[LB_PSI][i] - stMap[RL_PSI][i]) >= M_PI) 
                {   
                    double bl = stMap[LB_PSI][i] + 2 * M_PI * (stMap[LB_PSI][i] < 0);
                    double p = stMap[RL_PSI][i] + 2 * M_PI * (stMap[RL_PSI][i] < 0);
                    psi_interp = bl + (p - bl) * node_idx / raceline_index;          
                }
                else {
                    psi_interp = stMap[LB_PSI][i] + (stMap[RL_PSI][i] - stMap[LB_PSI][i]) * (node_idx+1) / raceline_index;
                }
                node.psi = normalizeAngle(psi_interp);
            }
            else if (node_idx == raceline_index) {
                psi_interp = stMap[RL_PSI][i];
                node.psi = psi_interp;
            }
            else {
                int remain = num_nodes - raceline_index - 1;
                double t = static_cast<double>(node_idx - raceline_index) / max(remain, 1);  // 0 ~ 1
                psi_interp = stMap[RL_PSI][i] + t * (stMap[RB_PSI][i] - stMap[RL_PSI][i]);
                node.psi = normalizeAngle(psi_interp);
            }
            // cout << i << "번째 레이어의" <<node_idx << "번째 노드의 psi는" << node.psi << endl;

            nodesPerLayer[i][node_idx] = node;
            ++node_idx;

        }

    }
    return {nodesPerLayer, raceline_index_array};
}

void calcOfflineCost(SplineMap& splineMap,
                   IVector& raceline_index_array,
                   YAML::Node& params) {
    if (splineMap.size() <= 0) {
        throw invalid_argument("SplineMap's Size is zero!!");
    }

    float w_raceline = params["cost"]["w_raceline"].as<float>();
    float w_raceline_sat = params["cost"]["w_raceline_sat"].as<float>();
    float w_length = params["cost"]["w_length"].as<float>();
    float w_curv_avg = params["cost"]["w_curv_avg"].as<float>();
    float w_curv_peak = params["cost"]["w_curv_peak"].as<float>();
    float lat_resolution = params["lattice"]["lat_resolution"].as<float>();

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
                cerr << "[WARNNING] Skipping spline: end_layer=" << end_layer 
              << " out of bounds (0.." << raceline_index_array.size()-1 << ")\n";
                continue;
            }

            if (spline.kappa.size() == 0)
            {
                // cerr << "[WARNNING] Skipping spline: empty curvature data\n";
                continue;
            }

            double abs_kappa = spline.kappa.array().abs().sum();
            double s_length = spline.el_lengths.sum();
            // cout << "s_length: " << s_length << endl;
            // average curvature
            offline_cost += w_curv_avg * pow(abs_kappa / float(spline.kappa.size()), 2) * s_length;
            // peak curvature
            double max_min = abs(spline.kappa.array().maxCoeff() - spline.kappa.array().minCoeff());
            offline_cost += w_curv_peak * pow(max_min, 2) * s_length;

            // path length
            offline_cost += w_length * s_length;

            // raceline cost

            double raceline_dist = abs(raceline_index_array[end_layer] - end_node) * lat_resolution;
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

void setInitialPose(DMap &stMap,
                    const NodeMap &nodesPerLayer,
                    const IVector &raceline_index_array,
                    YAML::Node &params) {
    // 현재 pos, heading 
    double dx, dy;
        
    Vector2d initial_pos(stMap[RL_X][1], stMap[RL_Y][1]);
    float vel_est = 0.0;

    float veh_width = params["vehicle"]["veh_width"].as<float>();
    double max_heading_offset = params["custom"]["max_heading_offset"].as<double>();
    float stepsize_approx = params["sampling"]["stepsize_approx"].as<float>();
    int initial_layer = params["planningtarget"]["initial_layer"].as<int>();
    
    // set start pos 
    if (!checkInsideBounds(stMap, initial_pos, veh_width)) {
        throw out_of_range("start pos is not in bounds");
    }
    
    IPair closest_idx;
    getClosestNodes(nodesPerLayer, closest_idx, initial_pos);
    int start_layer = closest_idx.first;
    int start_node = closest_idx.second;
    double start_heading = nodesPerLayer[closest_idx.first][closest_idx.second].psi;

    int end_layer = (closest_idx.first + initial_layer) % (nodesPerLayer.size() - 1);

    for (int layer_idx = start_layer; layer_idx < end_layer;++layer_idx) {
        if (layer_idx != start_layer) {
            start_node = raceline_index_array[layer_idx];
            start_heading = nodesPerLayer[layer_idx][start_node].psi;
        }

        int goal_layer = layer_idx + 1;
        int goal_node = raceline_index_array[goal_layer];
        double goal_heading = nodesPerLayer[goal_layer][goal_node].psi;
        double heading_diff = abs(start_heading - goal_heading);
        Vector2d start_pos(nodesPerLayer[layer_idx][start_node].x, nodesPerLayer[layer_idx][start_node].y);
        Vector2d end_pos(nodesPerLayer[goal_layer][goal_node].x, nodesPerLayer[goal_layer][goal_node].y);

        if (heading_diff > M_PI) {
            heading_diff = abs(2*M_PI - heading_diff);

        }
        if (heading_diff > max_heading_offset) {
            // cout << "heading_diff: " << heading_diff << ", " << max_heading_offset << endl;
            // cerr << "Heading mismatch between vehicle and track grid, check if vehicle oriented correctly!" << endl;
        }

        MatrixXd path(2, 2);
        path.block<1, 2>(0, 0) = start_pos.transpose();
        path.block<1, 2>(1, 0) = end_pos.transpose();

        auto result = calcSplines(path, start_heading, goal_heading); // coeffs_x, coeffs_y, kappa, el_lengths, cost
        auto [kappa, psi] = interpSplines(stMap, result->coeffs_x, result->coeffs_y, stepsize_approx, veh_width);

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

    string yaml_path = "config/offline_params.yaml";
    YAML::Node params = YAML::LoadFile(yaml_path);
    
    // unique_ptr<string> track = Load("include/driving_task.ini");
    // 자동 경로 설정
    string map_file_in = "maps/" + params["map_name"].as<string>() + ".csv";
    string map_file_out = "outputs/"+ params["map_name"].as<string>() + "_out.csv";

    // global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    // 결과: gtpl_map
    DMap gtMap = readDMapFromCSV(map_file_in);

    // 결과: gtpl_map에 삽입
    auto [rb_x, rb_y] = computeBoundRight(gtMap[POS_X], gtMap[POS_Y],
                                          gtMap[NORM_X], gtMap[NORM_Y],
                                          gtMap[WIDTH_R]);
    auto [lb_x, lb_y] = computeBoundLeft(gtMap[POS_X], gtMap[POS_Y],
                                         gtMap[NORM_X], gtMap[NORM_Y],
                                         gtMap[WIDTH_L]);
    auto [rl_x, rl_y] = computeRaceline(gtMap[POS_X], gtMap[POS_Y],
                                        gtMap[NORM_X], gtMap[NORM_Y],
                                        gtMap[NORM_L]);
    DVector rl_ds = computeDeltaS(gtMap[RL_S]);

    gtMap[RB_X] = rb_x;
    gtMap[RB_Y] = rb_y;
    gtMap[LB_X] = lb_x;
    gtMap[LB_Y] = lb_y;
    gtMap[RL_X] = rl_x;
    gtMap[RL_Y] = rl_y;
    gtMap[RL_dS] = rl_ds;

    // 결과: map_file_out 
    // [지민] gtpl_mp은 이미 전역변수라 삽입해줄 필요가 없음.
    writeDMapToCSV(map_file_out, gtMap);
    
    // layer 간격을 위한 raceline points sampling 
    IVector idx_sampling = samplePointsFromRaceline(gtMap[RL_KAPPA],
                                                    gtMap[RL_dS],
                                                    params);
    
    // cout << "idx size:" << idx_sampling.size() << endl;
    // 샘플링한대로 이제 gtpl_map 대신 stMap으로 바굼
    DMap stMap;
    for (const auto& [key, vec] : gtMap) {
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < vec.size()) {
                stMap[key].push_back(vec[idx]);
            }
        }
    }
    // writeDMapToCSV("inputs/stMap", stMap);
    // map_size(stMap); // (51, 3)

    // 결과: stMap에 delta_s열 추가 
    stMap[RL_dS] = computeDeltaS(stMap[RL_S]);

    // 결과: stMap에 열 추가
    stMap[RL_PSI] = calcHeading(stMap[RL_X], stMap[RL_Y]);
    // 여기서 계산되는 sampling된 bound_l, r은 node 생성 시에만 쓰인다. 
    stMap[LB_PSI] = calcHeading(stMap[LB_X], stMap[LB_Y]);
    stMap[RB_PSI] = calcHeading(stMap[RB_X], stMap[RB_Y]);  
    
    auto [nodesPerLayer, raceline_index_array] = genNode(stMap, params);
    
    // writeDMapToCSV("inputs/stMap.csv", stMap);

    auto [wayptGraph, splineMap] = genEdges(stMap, nodesPerLayer, raceline_index_array, params);
    
    // 결과: splineMap의 spline 구조체에 cost저장 
    calcOfflineCost(splineMap,
                   raceline_index_array,
                   params);
    
    // 결과: 초기경로 시각화
    setInitialPose(stMap,
                   nodesPerLayer,
                   raceline_index_array,
                   params);

    f_time = clock();

    // wayptGraph.printGraph();

    // visual process 
    cout << (double)(f_time - s_time) / CLOCKS_PER_SEC << "s 소요" << endl;
    
    // printSplineInfo(splineMap, nodesPerLayer);
    // 결과: 시각화
    visual(gtMap, stMap, nodesPerLayer, splineMap);

    return 0;
}