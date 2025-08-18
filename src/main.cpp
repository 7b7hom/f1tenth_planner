#include "spline.hpp"
// #include "graph_planner.hpp"

void visual(DMap &gtMap, DMap &stMap, const NodeMap &nodesPerLayer, SplineMap &splineMap);
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

bool checkInsideBounds(DMap &stMap, const Vector2d& pos, const float veh_width) {

    if (stMap.find(LB_X) == stMap.end() || 
    stMap.find(LB_Y) == stMap.end() ||
    stMap.find(RB_X) == stMap.end() || 
    stMap.find(RB_Y) == stMap.end()) {
    throw invalid_argument("Boundary keys are missing in stMap!");
}

    int n = stMap[LB_X].size();
    MatrixXd bound_l(n,2);
    MatrixXd bound_r(n,2);
    for (int i = 0; i < n; ++i) {
        bound_l(i, 0) = stMap[LB_X][i];
        bound_l(i, 1) = stMap[LB_Y][i];

        bound_r(i, 0) = stMap[RB_X][i];
        bound_r(i, 1) = stMap[RB_Y][i];
    }
    
    MatrixXd centerline = (bound_l + bound_r) / 2;

    // 가장 가까운 segment 인덱스 찾기
    int closest_idx = -1;
    double min_dist2 = numeric_limits<double>::max();
    for (int i = 0; i < centerline.rows() - 1; ++i) {
        // segment 중심 계산
        Vector2d mid = (centerline.row(i) + centerline.row(i + 1)) / 2.0;
        double dist2 = (mid - pos).squaredNorm();
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closest_idx = i;
        }
    }

    if (closest_idx < 0 || closest_idx >= bound_l.rows() - 1)
        return false; // 예외 처리

    // bound_l, bound_r, centerline 보간 (선형 보간 10개 지점)
    int interp_points = 10;
    MatrixXd bl_interp(interp_points, 2);
    MatrixXd br_interp(interp_points, 2);
    MatrixXd center_interp(interp_points, 2);

    for (int i = 0; i < interp_points; ++i) {
        double t = static_cast<double>(i) / (interp_points - 1);
        bl_interp.row(i) = (1 - t) * bound_l.row(closest_idx) + t * bound_l.row(closest_idx + 1);
        br_interp.row(i) = (1 - t) * bound_r.row(closest_idx) + t * bound_r.row(closest_idx + 1);
        center_interp.row(i) = (1 - t) * centerline.row(closest_idx) + t * centerline.row(closest_idx + 1);
    }

    // pos에 가장 가까운 center_interp 인덱스 찾기
    int nearest_idx = -1;
    double best_dist2 = numeric_limits<double>::max();
    for (int i = 0; i < interp_points; ++i) {
        double d2 = (center_interp.row(i) - pos.transpose()).squaredNorm();
        if (d2 < best_dist2) {
            best_dist2 = d2;
            nearest_idx = i;
        }
    }

    // bound 사이 거리 (제곱)
    double d_track2 = (bl_interp.row(nearest_idx) - br_interp.row(nearest_idx)).squaredNorm();

    // 차량에서 각 bound까지 거리 (제곱)
    double d_bl_2 = (bl_interp.row(nearest_idx) - pos.transpose()).squaredNorm();
    double d_br_2 = (br_interp.row(nearest_idx) - pos.transpose()).squaredNorm();

    double dist_to_left_bound = sqrt(d_bl_2);
    double dist_to_right_bound = sqrt(d_br_2);


    // cout << "-------here" << endl;
    // cout << dist_to_left_bound << endl;
    // cout << dist_to_right_bound << endl;
    // VEH_WIDTH 조건 확인
    if (dist_to_left_bound < veh_width || dist_to_right_bound < veh_width)
    {
        // throw invalid_argument("Spline point violates VEH_WIDTH constraints!");
        return false;
    }

    // bound 밖에 있는지 여부 확인
    bool within_bounds = !(d_bl_2 > d_track2 || d_br_2 > d_track2);
    return within_bounds;
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
    SplineHandler handler;

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

        auto result = handler.calcSplines(path, start_heading, goal_heading); // coeffs_x, coeffs_y, kappa, el_lengths, cost
        auto [points_xy, kappa] = handler.samplingSpline(result->coeffs_x, result->coeffs_y, params);

        if (kappa.size() == 0 || points_xy.size() == 0) {
            throw invalid_argument("Points or Kappa's Size is zero!! - interpSplines()");
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
    SplineHandler splineMap;

    auto wayptGraph = splineMap.genEdges(stMap, nodesPerLayer, raceline_index_array, params);
    cout << "Initial generated splines: ";
    wayptGraph.printGraph();
    splineMap.writeSplineMapToCSV("outputs/splineMap.csv");

    float veh_turn = params["vehicle"]["veh_turn"].as<float>();
    float min_vel_race = params["lattice"]["min_vel_race"].as<float>();
    float max_lateral_accel = params["lattice"]["max_lateral_accel"].as<float>();
    float veh_width = params["vehicle"]["veh_width"].as<float>();
    int rmv_cnt = 0;

    for (size_t layer_idx = 0; layer_idx < nodesPerLayer.size();++layer_idx) {
      int srcLayerIdx = layer_idx;
      for (size_t node_idx = 0; node_idx < nodesPerLayer[srcLayerIdx].size(); ++node_idx) {
        IPairVector childNodes;
        IPair start = make_pair(layer_idx ,node_idx);
        wayptGraph.getChildNodes(start, childNodes);

        // 연결된 노드와 loop
        for (auto& end : childNodes) {
            MatrixXd& coeffs_x = splineMap.at(start, end).coeffs_x;
            MatrixXd& coeffs_y = splineMap.at(start, end).coeffs_y;
            // spline 위의 점들을 샘플링(no_interp_points개수만큼)
            auto [points_xy, kappa] = splineMap.samplingSpline(coeffs_x, coeffs_y, params);
            // 점들을 기준으로 pruneEdges에 가서 1.곡률 2.트랙내 여부 에 따라 remove를 한다.
            if (kappa.size() == 0 || points_xy.size() == 0) {
                cerr << "Invalid spline sampling" << endl;
                continue;
            }
            // 해당 spline위에서 샘플링한 점이 track을 벗어나면 pruneEdges()에 갈 수 있도록.
            splineMap.at(start, end).kappa = kappa;
            splineMap.at(start, end).points_xy = points_xy;

            if (!splineMap.at(start, end).raceline) {

                int layer_idx = start.first;
                double vel_rl = stMap[RL_VX][layer_idx] * min_vel_race;
                double min_turn = pow(vel_rl, 2) / max_lateral_accel;

                bool toRemove = false;

                for (int j = 0; j < kappa.size(); ++j) {

                    double kappa_val = abs(kappa(j));
                    
                    if ((kappa_val > 1.0 / veh_turn || kappa_val > 1.0 / min_turn))
                    {
                        toRemove = true;
                        break;
                    }
                }
                if (toRemove) {
                    wayptGraph.removeEdge(start, end, &splineMap.getSplineMap(), static_cast<int>(nodesPerLayer.size()));
                    rmv_cnt++;
                    }
                }
            }
        }
    }
    cout << "Number of splines deleted due to curvature conditions: " << rmv_cnt << endl;

    splineMap.pruneEdge(wayptGraph, nodesPerLayer);

    cout << "The number of splines generated finally: ";
    wayptGraph.printGraph();

    // 결과: splineMap의 spline 구조체에 cost저장 
    splineMap.calcOfflineCost(raceline_index_array, params);
    f_time = clock();
    
    // 결과: 초기경로 시각화
    setInitialPose(stMap,
                   nodesPerLayer,
                   raceline_index_array,
                   params);

    // wayptGraph.printGraph();

    // visual process 
    cout << (double)(f_time - s_time) / CLOCKS_PER_SEC << "s 소요" << endl;

    splineMap.readSplineMapFromCSV("splineMap.csv");
    // printSplineInfo(splineMap, nodesPerLayer);
    // 결과: 시각화
    visual(gtMap, stMap, nodesPerLayer, splineMap.getSplineMap());

    return 0;
}