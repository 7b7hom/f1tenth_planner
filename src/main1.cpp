    string map_file_out = "inputs/gtpl_levine_out.csv";
#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// Dvector를 Map 구조로 추가(연산)
void addDVectorToMap(DMap &map,
                     string attr,
                     const IVector *idx_array = nullptr) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();
    } 
    else {
        len = idx_array->size();
    }
    // cout << "attr: "<< attr << " / len:" << len << endl;

    DVector x_out(len), y_out(len);
    string x_label = "x_" + attr;
    string y_label = "y_" + attr;
    
    if (!attr.compare("bound_r")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__width_right][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__width_right][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("bound_l")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] - map[__x_normvec][i] * map[__width_left][i];
            y_out[i] = map[__y_ref][i] - map[__y_normvec][i] * map[__width_left][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("raceline")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__alpha][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__alpha][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    // i번째와 i-1번째 point의 delta_s 계산 
    // delta_s[0] = 0 
    else if (!attr.compare("delta_s")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len - 1; ++i) {
            x_out[i] = map[__s_racetraj][i+1] - map[__s_racetraj][i]; // 마지막 원소는 0
        }
        map[attr] = x_out; 
    }

    // map_size(map);
}

void samplePointsFromRaceline(const DVector& kappa,
                              const DVector& dist,
                              double d_curve,
                              double d_straight,
                              double curve_th,
                              IVector& idx_array) {

    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    for (size_t i = 0; i < n; ++i) {

        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist;
        }

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

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

void calcHeading(DVector &x_raceline,
                 DVector &y_raceline,
                 DVector &psi) {

    size_t N = x_raceline.size();
    psi.resize(N);

    // 닫힌 회로 가정. 예외 처리 필요
    double dx, dy;
    for (size_t i = 0; i < N; ++i) {
        
        if (i != N -1) {
            dx = x_raceline[i+1] - x_raceline[i];
            dy = y_raceline[i+1] - y_raceline[i];
        } else {
            dx = x_raceline[0] - x_raceline[N - 1];
            dy = y_raceline[0] - y_raceline[N - 1];
        } 
    psi[i] = atan2(dy, dx) - M_PI_2;
        
    normalizeAngle(psi[i]);

    }
    // cout << i<< ": " << psi[i] << endl;
    // cout << psi.size() << endl;

}

void genNode(NodeMap& nodesPerLayer,
            const double veh_width,
            float lat_resolution) {
    
    const size_t N = sampling_map[__alpha].size();
    IVector raceline_index_array;
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        node.layer_idx = i; 
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);
        
        cout << "layer 길이" << (sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width/2)<< endl;
        cout << "layer 내에서 raceline index:" << raceline_index << endl;
        cout << "-----" << endl;

        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);    // 기준선에서의 위치
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);  // 기준선에서 수직한 노멀 벡터 따라 노드 배치
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지
        int node_idx = 0;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;  // num_nodes : 좌우 총 가능한 노드 수
        nodesPerLayer[i].resize(num_nodes); 

        // cout << i << "번째 layer의 node 개수는 " << num_nodes << endl;
        // node별 loop 
        for (double alpha = start_alpha; alpha <= sampling_map[__width_right][i] - veh_width / 2 ; alpha+=lat_resolution) {
            // node_alphas.push_back(alpha);
            // node의 좌표 계산.
            node_pos = ref_xy + alpha * norm_vec;
            // node의 layer내의 인덱스 계산.
            node.node_idx = node_idx;
            node.x = node_pos.x();
            node.y = node_pos.y();      
            node.raceline = (node_idx == raceline_index);
            
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
        // 곡률 계산 (헤딩 변화량 / 거리)
        for (int j = 1; j < num_nodes - 1; ++j) {
            const Node& prev = nodesPerLayer[i][j - 1];
            const Node& next = nodesPerLayer[i][j + 1];

            double dpsi = normalizeAngle(next.psi - prev.psi);
            double ds = std::hypot(next.x - prev.x, next.y - prev.y);
            double kappa = (ds > 1e-6) ? dpsi / ds : 0.0;
            nodesPerLayer[i][j].kappa = kappa;
        }

        // 양 끝은 0으로 처리
        nodesPerLayer[i][0].kappa = 0.0;
        nodesPerLayer[i][num_nodes - 1].kappa = 0.0;

        // cout << i << "번째 Layer의" << endl;
        // for (size_t i =0; i < node_pos.size(); ++i) {        
        //     cout << i << "번째 Node" << endl;
        //     cout << node_pos[i] << endl;
        // }

    }
}

void calcSplines() {

}

void genEdges(NodeMap &nodesPerLayer, 
              IVector &raceline_index_array,
              Vector2d &node_pos,
              Graph &edgeList,
              const double lat_offset,
              const double lat_resolution,
              const double const_thr,
              const double min_vel = 0.0,
              bool closed = true) {

    if (lat_offset <= 0.0) {
        throw invalid_argument("Too small lateral offset!");
    }

    vector<Vector2d> raceline_cl;
    Vector2d start_point;

    if (closed) {
        raceline_cl.push_back(node_pos.x(), node_pos.y());
    }
    // raceline spline 먼저 계산
    auto raceline_result = calcSplines(raceline_cl);

    const MatrixXd& x_coeff_r = raceline_result -> coeffs_x;
    const MatrixXd& y_coeff_r = raceline_result -> coeffs_y;

    // 레이어 별 loop
    for (size_t layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) {
        
        int start_layer = layerIdx;
        int end_layer = layerIdx + 1;

        // 마지막 layer의 경우 0번째 layer와 연결시킬 수 있도록 end_layer 조정 
        if (end_layer >= nodesPerLayer.size())
            end_layer -= nodesPerLayer.size();
        else 
            break;

        int start_raceline_Idx = raceline_index_array[start_layer];
        int end_raceline_Idx = raceline_index_array[end_layer];
        
        // start layer 내 노드별 loop
        for (size_t startIdx = 0; startIdx <= nodesPerLayer[startLayer].size(); ++startIdx) {
            // 기준 노드
            Node &startNode = nodesPerLayer[startLayer][startIdx];

            int refDestIdx = end_raceline_Idx + startIdx - start_raceline_Idx; // 기준 end노드 
            
            refDestIdx = clamp(refDestIdx, 0, nodesPerLayer[end_layer].size() - 1);
            Node &srcEndNode = nodesPerLayer[end_layer][refDestIdx];
            Vector2d d_start(startNode.x, startNode.y);
            Vector2d d_end(srcEndNode.x, srcEndNode.y);

            // spline 연결할 노드 선정 기준 : lat_steps
            double dist = (d_end - d_start).norm();
            // genNode에서 kappa 계산한거 토대로(+기능 추가 필요)
            double factor = (startNode.kappa > curve_thr) ? 2.0 : 1.0;  // 커브에서 더 많이 연결(추월 경로를 위하여)
            int lat_steps = round(factor * dist * lat_offset / lat_resolution); // srcEndNode 기준 2*lat_steps + 1개의 노드와 연결한다.

            // startNode와 lat_steps 기준 해당되는 노드들 spline 연결 
            for (int destIdx = max(0, refDestIdx - lat_steps); 
                desNode <= min(nodesPerLayer[end_layer].size() -1, refDestIdx + lat_steps); ++destIdx) {
                    Node &endNode = nodesPerLyaer[end_layer][destIdx];

                    MatrixXd x_coeffs, y_coeffs;

                    if (srcNode.raceline && endNode.raceline ) {
                        x_coeffs = x_coeff_r.row(start_layer);
                        y_coeffs = y_coeff_r.row(start_layer);
                    }
                    else {
                    vector<Vector2d> path = { Vector2d(start_node.x, start_node.y), Vector2d(end_node.x, end_node.y) };
                    auto result = calcSplines(path, nullptr, start_node.psi, end_node.psi);
                    x_coeff = result->coeffs_x;
                    y_coeff = result->coeffs_y;
                    }
                    // graph에 넣는 과정 
                    IPair n1 = make_pair(start_layer, startIdx);
                    IPair n2 = make_pair(end_layer, destIdx);
                    edgeList.addEdge(n1, n2);

                }
        }
    }
}

int main() {
    IVector idx_sampling;
    Offline_Params params;

    string map_file_in = "inputs/gtpl_levine.csv";
    string map_file_out = "inputs/gtpl_levine_out.csv";

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

    genNode(nodesPerLayer,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    // sampling points' info 
    // writeDMapToCSV("inputs/sampling_map.csv", sampling_map);
    // visual process 
    visual(nodesPerLayer);

    #if 0
    // Graph sample code 
    Graph directedGraph;
    IPair t1 = make_pair(0, 0);
    IPair t2 = make_pair(0, 1);
    IPair t3 = make_pair(1, 0);
    IPair t4 = make_pair(0, 2);
    IPair t5 = make_pair(1, 2);
    IPair t6 = make_pair(1, 5);
    IPair t7 = make_pair(1, 6);
    // spline 생성 후 edge로 집어넣음.
    directedGraph.addEdge(t1, t3);// push_back이라서 sorting은 되지 않음. 
    directedGraph.addEdge(t1, t5);
    directedGraph.addEdge(t1, t6);
    directedGraph.addEdge(t2, t5);
    directedGraph.addEdge(t4, t7);
    // 실제 로직은 node idx가 작은 순서대로 그래프가 그러질 예정이라 괜찮을 듯.
    // grpah 전체 print 
    cout << "---처음 Graph---" << endl;
    directedGraph.printGraph();

    IPairVector child1;
    // t1 노드의 뒤로 연결된(child) node들을 뽑아온다.
    // directedGraph.getChildIdx(t1, child1);

    // for (size_t i =0; i < child1.size(); ++i) {
    //     cout << child1[i].first << ", " << child1[i].second << "/ ";
    // }
    
    // Error Index(child가 없는 경우 runtime_error)
    // directedGraph.getChildIdx(t5, child1);
    
    // (0, n)이라는 임의의 노드 n이 t5(1, 2)을 들고 있는 경우 해당 list에서 t5 삭제
    // 이후 child 노드가 그 뒤로 연결된 spline이 없는 경우 parent의 adjList에서 child 노드를 삭제하기 위하여 필요함.
    IPairVector parent; 
    // t5 노드를 들고 있는 노드가 있는지 1. 찾고 2. parent로 반환함.
    directedGraph.getParentNode(t5, parent); 
    
    cout << "---0번째 Layer의 node 중에서 1번째 Layer의 3번째 노드와 엣지로 연결되어 있는 노드의 idx---" << endl;
    for (size_t i = 0; i < parent.size(); ++i) {
        cout << parent[i].first << ", " << parent[i].second << endl;
 
    }
    directedGraph.removeEdge(t5, parent);
    cout << "---위의 엣지를 제거한 후 graph 상태---" << endl;
    // 결과 확인용 
    directedGraph.printGraph();
    #endif

    return 0;
}