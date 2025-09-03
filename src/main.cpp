#include "graph_planner.hpp"
#include "config_modena.h"

#include <chrono>
using clock_type = std::chrono::steady_clock;

// 전체 경로 계획 파이프라인을 실행하는 함수
int main() {
    Offline_Params params;
    auto T0 = clock_type::now();

    string map_file_in = "/home/uiiiqns/f1tenth_planner/inputs/traj_ltpl_cl_modena.csv";
    string map_file_out = "/home/uiiiqns/f1tenth_planner/inputs/traj_ltpl_cl_modena_out.csv";

    // 1. 트랙 데이터 로드 및 전처리
    auto t = clock_type::now();
    DMap gtpl_map = readDMapFromCSV(map_file_in);
    auto t1 = clock_type::now();
    std::cout << "[T] read csv: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n";

    // CSV로부터 반드시 있어야 하는 열들(스플라인/경계 계산 전에)
    assert_has_keys(gtpl_map,
        {__x_ref, __y_ref, __x_normvec, __y_normvec, __width_right, __width_left, __alpha,
        __kappa, __s_racetraj}, "gtpl_map after readDMapFromCSV");

    // 1-1. 계산
    auto [x_right, y_right] = computeBoundRight(gtpl_map);
    auto [x_left, y_left] = computeBoundLeft(gtpl_map);
    auto [x_race, y_race] = computeBoundRace(gtpl_map);
    DVector ds_race = computeDeltaS(gtpl_map[__s_racetraj], true);

    // 1-2. DMap에 등록
    gtpl_map[__x_bound_r] = move(x_right);
    gtpl_map[__y_bound_r] = move(y_right);
    gtpl_map[__x_bound_l] = move(x_left);
    gtpl_map[__y_bound_l] = move(y_left);
    gtpl_map[__x_raceline] = move(x_race);
    gtpl_map[__y_raceline] = move(y_race);
    gtpl_map[__delta_s] = move(ds_race);
 
    writeDMapToCSV(map_file_out, gtpl_map, ',');

    t = clock_type::now();

    // 2. 레이어 샘플링
    IVector idx_sampling = samplePointsFromRaceline(gtpl_map[__kappa], gtpl_map[__delta_s],
                             params.d_curve, params.d_straight,
                             params.curve_thr);

    DMap sampled_map;
    
    // 샘플링된 인덱스를 사용하여 gtpl_map에서 데이터를 복사하여 sampling_map 채우기
    for(const auto& kv : gtpl_map){
        const auto& key = kv.first; 
        const auto& vec = kv.second;
        DVector col; col.reserve(idx_sampling.size());
        for(int idx : idx_sampling) if (0 <= idx && idx < (int)vec.size()) col.push_back(vec[idx]);
        sampled_map[key] = move(col);
    }
    // addDVectorToMap(sampling_map, "delta_s", &idx_sampling); 

    t1 = clock_type::now();
    std::cout << "[T] sampling + copy: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
            << " ms\n";

    // 샘플링 끝나고 헤딩 계산 전
    assert_has_keys(sampled_map, {__x_raceline, __y_raceline, __x_bound_l, __y_bound_l, __x_bound_r, __y_bound_r},
                    "sampled_map before heading");
    
    t = clock_type::now();
    sampled_map[__psi] = calcHeading(sampled_map[__x_raceline], sampled_map[__y_raceline]);
    sampled_map[__psi_bound_l] = calcHeading(sampled_map[__x_bound_l], sampled_map[__y_bound_l]);
    sampled_map[__psi_bound_r] = calcHeading(sampled_map[__x_bound_r], sampled_map[__y_bound_r]);
    t1 = clock_type::now();
    std::cout << "[T] headings: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n";

    // 3. 노드 그리드 생성
    t = clock_type::now();
    auto layers = computeLayerParams(sampled_map, params.veh_width, params.lat_resolution);
    NodeMap nodes = buildNodeGrid(layers, sampled_map, params.lat_resolution);
    fillNodeHeadings(nodes, layers, sampled_map);
    t1 = clock_type::now();
    std::cout << "[T] nodes grid+heading: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n";

    //NodeMap nodesPerLayer = genNode(sampled_map, params.veh_width, params.lat_resolution);

    // 4. 스플라인 생성 및 유효성 검사, 최종 그래프 구축
    t = clock_type::now();
    Graph graph = makeEmptyGraph(nodes);
    addRacelineEdges(graph, nodes, params, sampled_map);
    addCandidateEdges(graph, nodes, params, sampled_map);
    t1 = clock_type::now();
    std::cout << "[T] edges (spline+check): "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n";

    //Graph graph = genEdges(nodesPerLayer, params, sampled_map);
    graph.printGraph(); 

    // 5. prune
    t = clock_type::now();
    graph = prune_graph(move(graph), (int)nodes.size(), true); 
    t1 = clock_type::now();
    std::cout << "[T] prune: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n";

    // 6. calc cost
    t = clock_type::now();
    graph = gen_offline_cost(move(graph), params, nodes);
    t1 = clock_type::now();
    std::cout << "[T] offline cost: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t).count()
              << " ms\n"; 

    cout << "\n--- 최종 생성된 그래프 (유효한 스플라인 엣지 포함) ---" << endl;
    graph.printGraph();

    // (7. set_startpos())

    Vector2d pos_est(143.45, -130.90); 
    double heading_est = -2.21; // rad

    // set_startpos 함수를 호출하여 가장 가까운 노드를 찾고 엣지를 생성합니다.
    StartPosResult spr = set_startpos(pos_est, heading_est, params, nodes,
                sampled_map, 20.0 * M_PI / 180.0);
                
    if(spr.out_of_track){
        cout << "Out of track or bad heading -> start not fixed, no edge created.\n";
    }else{
        graph = ensure_edge_between(move(graph), nodes, params, spr.start_key, spr.next_idx, sampled_map);
        int l0 = get<0>(spr.start_key), i0 = get<1>(spr.start_key);
        int l1 = (l0 + 1) % (int)nodes.size();
        cout << "[FIXED] (" << l0 << "," << i0 << ")->(" << l1 << "," << spr.next_idx << ")\n";

        const auto& adj = graph.getAdjLists();
        bool has_edge = false;
        if(auto it = adj.find(spr.start_key); it != adj.end()){
            has_edge = (it->second.count(spr.next_idx) > 0);
        }
        cout << "Edge created?" << (has_edge ? "YES" : "NO") << "\n";
    }

    auto T1 = clock_type::now();
    std::cout << "[T] TOTAL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()
              << " ms\n";

    // 8. 결과 시각화 (Plotting)
    //visual(gtpl_map, sampled_map, nodes, graph, params, pos_est, heading_est, spr.start_key, spr.next_idx);
}
