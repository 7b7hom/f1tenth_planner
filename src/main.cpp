#include "graph_planner.hpp"
#include "config_modena.h"

extern DMap gtpl_map;
extern DMap sampling_map;

void readDMapFromCSV(const std::string& pathname, DMap& map);
void writeDMapToCSV(const std::string& pathname, DMap& map, char delimiter);
void addDVectorToMap(DMap &map, std::string attr, const IVector *idx_array);
void samplePointsFromRaceline(const DVector& kappa, const DVector& dist,
                              double d_curve, double d_straight, double curve_th, IVector& idx_array);
void calcHeading(DVector &x, DVector &y, DVector &psi);

void genNode(NodeMap& nodesPerLayer, const double veh_width, float lat_resolution);
void genEdges(Graph& graph, const NodeMap& nodesPerLayer, const Offline_Params params);
void prune_graph(Graph& graph, int num_layers, bool closed);
void gen_offline_cost(Graph& graph, const Offline_Params& params, const NodeMap& nodesPerLayer);

void set_startpos(const Eigen::Vector2d& pos_est, double heading_est, const Offline_Params& params,
                  const NodeMap& nodesPerLayer, Graph& graph, double max_heading_offset_rad,
                  ITuple& out_start_key, int& out_next_idx, bool& out_of_track);

void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params,
            const Vector2d& pos_est, double heading_est,
            const ITuple& start_key, int next_idx);

// 전체 경로 계획 파이프라인을 실행하는 함수
int main() {
    Offline_Params params;

    std::string map_file_in = "/home/subin/subin/AID/planner_project/f1tenth_planner/inputs/traj_ltpl_cl_modena.csv";
    std::string map_file_out = "/home/subin/subin/AID/planner_project/f1tenth_planner/inputs/traj_ltpl_cl_modena_out.csv";

    // 1. 트랙 데이터 로드 및 전처리
    readDMapFromCSV(map_file_in, gtpl_map); // gen_spline.cpp의 전역 gtpl_map에 로드
    addDVectorToMap(gtpl_map, "bound_r", nullptr);
    addDVectorToMap(gtpl_map, "bound_l", nullptr);
    addDVectorToMap(gtpl_map, "raceline", nullptr);
    addDVectorToMap(gtpl_map, "delta_s", nullptr);
    writeDMapToCSV(map_file_out, gtpl_map, ',');

    // 2. 레이어 샘플링
    IVector idx_sampling;
    samplePointsFromRaceline(gtpl_map[__kappa], gtpl_map[__delta_s],
                             params.d_curve, params.d_straight,
                             params.curve_thr, idx_sampling);
    
    // 샘플링된 인덱스를 사용하여 gtpl_map에서 데이터를 복사하여 sampling_map 채우기
    for (const auto& [key, vec] : gtpl_map) {
        sampling_map[key].reserve(idx_sampling.size()); 
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < (int)vec.size()) { 
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }
    addDVectorToMap(sampling_map, "delta_s", &idx_sampling); 
    
    calcHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);
    calcHeading(sampling_map[__x_bound_l], sampling_map[__y_bound_l], sampling_map[__psi_bound_l]);
    calcHeading(sampling_map[__x_bound_r], sampling_map[__y_bound_r], sampling_map[__psi_bound_r]);

    // 3. 노드 그리드 생성
    NodeMap nodesPerLayer;
    genNode(nodesPerLayer, params.veh_width, params.lat_resolution);

    // Offline_Params에 총 레이어 수 설정 (gen_offline_cost와 prune_graph에서 사용)
    Offline_Params mutable_params = params; 
    //mutable_params.num_layers = nodesPerLayer.size();

    // 4. 스플라인 생성 및 유효성 검사, 최종 그래프 구축
    Graph directedGraph(nodesPerLayer.size()); 
    genEdges(directedGraph, nodesPerLayer, mutable_params); 

    directedGraph.printGraph();

    // 5. prune
    prune_graph(directedGraph, nodesPerLayer.size(), true); 

    // 6. calc cost
    gen_offline_cost(directedGraph, mutable_params, nodesPerLayer); 

    cout << "\n--- 최종 생성된 그래프 (유효한 스플라인 엣지 포함) ---" << endl;
    directedGraph.printGraph();

    // (7. set_startpos())
    ITuple start_key;
    int next_idx;
    bool out_of_track;

    std::cout << "First raceline heading: " << sampling_map[__psi][0] << " rad" << std::endl;
    Vector2d pos_est(143.45, -130.90); 
    double heading_est = -2.21; // rad

    // set_startpos 함수를 호출하여 가장 가까운 노드를 찾고 엣지를 생성합니다.
    set_startpos(pos_est, heading_est, mutable_params, nodesPerLayer,
                directedGraph, 20.0 * M_PI / 180.0,
                start_key, next_idx, out_of_track);

    if (out_of_track) {
        std::cout << "Out of track or bad heading -> start not fixed, no edge created.\n";
    } else {
        int l0 = get<0>(start_key), i0 = get<1>(start_key);
        int l1 = (l0 + 1) % nodesPerLayer.size();
        std::cout << "[FIXED] ("<< l0 <<","<< i0 <<")->("<< l1 <<","<< next_idx <<")\n";

        const auto& adj = directedGraph.getAdjLists();
        bool has_edge = false;
        if (auto it = adj.find(start_key); it != adj.end()) {
            has_edge = (it->second.count(next_idx) > 0);
        }
        std::cout << "Edge created? " << (has_edge ? "YES" : "NO") << "\n";
    }

    // 8. 결과 시각화 (Plotting)
    visual(nodesPerLayer, directedGraph, params, pos_est, heading_est, start_key, next_idx);
}
