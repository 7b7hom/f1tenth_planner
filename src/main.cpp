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

void computeCurvature(NodeMap& nodesPerLayer) {
    const int num_layers = nodesPerLayer.size();

    // safe node access: 없는 j 인덱스는 가장 가까운 노드로 fallback
    auto safeGetNode = [](const std::vector<Node>& layer, int j) -> const Node* {
        if (layer.empty()) return nullptr;
        if (j < 0) return &layer.front();
        if (j < static_cast<int>(layer.size())) return &layer[j];
        return &layer.back();  // 가장 오른쪽 노드로 fallback
    };

    for (int i = 0; i < num_layers; ++i) {
        int num_nodes_i = nodesPerLayer[i].size();

        for (int j = 0; j < num_nodes_i; ++j) {
            const Node* prev = nullptr;
            const Node* next = nullptr;

            if (i > 0 && i < num_layers - 1) {
                prev = safeGetNode(nodesPerLayer[i - 1], j);
                next = safeGetNode(nodesPerLayer[i + 1], j);
            } else if (i == 0 && num_layers > 1) {
                prev = safeGetNode(nodesPerLayer[i], j);
                next = safeGetNode(nodesPerLayer[i + 1], j);
            } else if (i == num_layers - 1 && num_layers > 1) {
                prev = safeGetNode(nodesPerLayer[i - 1], j);
                next = safeGetNode(nodesPerLayer[i], j);
            }

            double dpsi = 0.0, ds = 0.0;
            if (prev && next) {
                dpsi = normalizeAngle(next->psi - prev->psi);
                ds = std::hypot(next->x - prev->x, next->y - prev->y);
            }

            double kappa = (ds > 1e-6) ? dpsi / ds : 0.0;
            nodesPerLayer[i][j].kappa = kappa;
        }
    }
}


void genNode(NodeMap& nodesPerLayer,
            IVector &raceline_index_array,
            const double veh_width,
            float lat_resolution) {
    
    const int N = sampling_map[__alpha].size();
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        node.layer_idx = i; 
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);
        
        // cout << "layer 길이" << (sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width/2)<< endl;
        // cout << "layer 내에서 raceline index:" << raceline_index << endl;
        // cout << "-----" << endl;

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

        
    }
    // 곡률 계산 (헤딩 변화량 / 거리)
    computeCurvature(nodesPerLayer);
        // cout << i << "번째 Layer의" << endl;
        // for (size_t i =0; i < node_pos.size(); ++i) {        
        //     cout << i << "번째 Node" << endl;
        //     cout << node_pos[i] << endl;
        // }

}

unique_ptr<SplineResult> calcSplines(
    const vector<Vector2d> &path,
    double psi_s = NAN,
    double psi_e = NAN) {

    vector<Vector2d> closed_path = path;
    bool closed = (path.front() - path.back()).norm() < 1e-6;

    // if (closed) closed_path.push_back(path[0]);

    int no_splines = closed_path.size() - 1;
    VectorXd ds = VectorXd::Ones(no_splines);

    MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
    VectorXd b_x = VectorXd::Zero(no_splines * 4);
    VectorXd b_y = VectorXd::Zero(no_splines * 4);

    for (int j = 0; j < no_splines; ++j) {
        int row0 = j * 4;
        M(row0, row0) = 1;
        b_x(row0) = closed_path[j].x();
        b_y(row0) = closed_path[j].y();

        RowVector4d T;
        double d = ds(j);
        T << 1, d, d * d, d * d * d;

        for (int k = 0; k < 4; ++k)
            M(row0 + 1, row0 + k) = T(k);

        b_x(row0 + 1) = closed_path[j + 1].x();
        b_y(row0 + 1) = closed_path[j + 1].y();
    }

    for (int j = 0; j < no_splines - 1; ++j) {
        int row = 2 * no_splines + j;
        double d = ds(j);
        RowVector4d T;
        T << 0, 1, 2 * d, 3 * d * d;
        for (int k = 0; k < 4; ++k) {
            M(row, 4 * j + k) = T(k);
            M(row, 4 * (j + 1) + k) = -T(k);
        }
    }

    for (int j = 0; j < no_splines - 1; ++j) {
        int row = 3 * no_splines - 1 + j;
        double d = ds(j);
        RowVector4d T;
        T << 0, 0, 2, 6 * d;
        for (int k = 0; k < 4; ++k) {
            M(row, 4 * j + k) = T(k);
            M(row, 4 * (j + 1) + k) = -T(k);
        }
    }

    if (!isnan(psi_s)) {
        M.row(4 * no_splines - 2).setZero();
        M(4 * no_splines - 2, 1) = 1;
        b_x(4 * no_splines - 2) = cos(psi_s);
        b_y(4 * no_splines - 2) = sin(psi_s);
    }

    if (!isnan(psi_e)) {
        double d = ds(no_splines - 1);
        M.row(4 * no_splines - 1).setZero();
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 1) = 1;
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 2) = 2 * d;
        M(4 * no_splines - 1, 4 * (no_splines - 1) + 3) = 3 * d * d;
        b_x(4 * no_splines - 1) = cos(psi_e);
        b_y(4 * no_splines - 1) = sin(psi_e);
    }

    VectorXd x_les = M.colPivHouseholderQr().solve(b_x);
    VectorXd y_les = M.colPivHouseholderQr().solve(b_y);

    MatrixXd coeffs_x = Map<MatrixXd>(x_les.data(), 4, no_splines).transpose();
    MatrixXd coeffs_y = Map<MatrixXd>(y_les.data(), 4, no_splines).transpose();

    MatrixXd normvec(no_splines, 2);
    for (int i = 0; i < no_splines; ++i) {
        double dx = coeffs_x(i, 1);
        double dy = coeffs_y(i, 1);
        normvec(i, 0) = -dy;
        normvec(i, 1) = dx;
    }

    VectorXd norms = normvec.rowwise().norm();
    MatrixXd normvec_normalized(no_splines, 2);
    for (int i = 0; i < no_splines; ++i)
        normvec_normalized.row(i) = normvec.row(i) / norms(i);

    return make_unique<SplineResult>(SplineResult{coeffs_x, coeffs_y, M, normvec_normalized});
}

void genEdges(NodeMap &nodesPerLayer, 
              IVector &raceline_index_array,
              Graph &edgeList,
              SplineMap &splineMap,
              const double lat_offset,
              const double lat_resolution,
              const double curve_thr,
              const double min_vel = 0.0,
              bool closed = true) {
    
    if (lat_offset <= 0.0) {
        throw invalid_argument("Too small lateral offset!");
    }

    // cout << nodesPerLayer.size() << endl; 출력: 51
    // 레이어 별 loop
    for (int layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) {
        
        int start_layer = layerIdx;
        int end_layer = layerIdx + 1;

        // cout << "end_layer:" << end_layer << endl;
        // cout << "nodesPerLayer.size()" << nodesPerLayer.size() << endl;

        // 마지막 layer의 경우 0번째 layer와 연결시킬 수 있도록 end_layer 조정 
        if (end_layer >= nodesPerLayer.size()) {
            end_layer -= nodesPerLayer.size();
        }

        int start_raceline_Idx = raceline_index_array[start_layer];
        int end_raceline_Idx = raceline_index_array[end_layer];

        // start layer 내 노드별 loop
        for (size_t startIdx = 0; startIdx <= nodesPerLayer[start_layer].size(); ++startIdx) {
            // 기준 노드
            Node &startNode = nodesPerLayer[start_layer][startIdx];

            // int refDestIdx = end_raceline_Idx + startIdx - start_raceline_Idx; // 기준 end노드 
            
            // refDestIdx = clamp(refDestIdx, 0, static_cast<int>(nodesPerLayer[end_layer].size() - 1));
            int refDestIdx = startIdx;
            Node &srcEndNode = nodesPerLayer[end_layer][refDestIdx];
            Vector2d d_start(startNode.x, startNode.y);
            Vector2d d_end(srcEndNode.x, srcEndNode.y);

            // spline 연결할 노드 선정 기준 : lat_steps
            double dist = (d_end - d_start).norm();
            // genNode에서 kappa 계산한거 토대로(+기능 추가 완료)
            double factor = (startNode.kappa > curve_thr) ? 2.0 : 1.0;  // 커브에서 더 많이 연결(추월 경로를 위하여)
            int lat_steps = round(factor * dist * lat_offset / lat_resolution); // srcEndNode 기준 2*lat_steps + 1개의 노드와 연결한다.

            // startNode와 lat_steps 기준 해당되는 노드들 spline 연결 
            for (int destIdx = max(0, refDestIdx - lat_steps); 
                destIdx <= min(static_cast<int>(nodesPerLayer[end_layer].size() -1), refDestIdx + lat_steps); ++destIdx) {
                    Node &endNode = nodesPerLayer[end_layer][destIdx];

                    vector<Vector2d> path = { Vector2d(startNode.x, startNode.y), Vector2d(endNode.x, endNode.y) };
                    auto result = calcSplines(path, startNode.psi, endNode.psi);

                    IPair startKey = make_pair(start_layer, startIdx);
                    IPair endKey = make_pair(end_layer, destIdx);
                    EdgeKey srcKey= make_pair(startKey, endKey);
                    splineMap[srcKey] = *result;

                    // graph에 넣는 과정 
                    edgeList.addEdge(startKey, endKey);

                    // cout << "startKey:" << startKey.first << ", " << startKey.second << " -> ";
                    // cout << "endKey:" << endKey.first << ", " << endKey.second << endl;
                    
                }
        }
    }
}

void printSplineMapVerbose(const SplineMap& splineMap, const NodeMap& nodesPerLayer) {
    cout << "\n=== 📌 SplineMap: Coefficients with Start/End Node Info ===\n";

    for (const auto& [edgeKey, spline] : splineMap) {
        const IPair& startKey = edgeKey.first;
        const IPair& endKey = edgeKey.second;

        const Node& startNode = nodesPerLayer[startKey.first][startKey.second];
        const Node& endNode = nodesPerLayer[endKey.first][endKey.second];

        cout << "\n▶ (" << startKey.first << ", " << startKey.second << ") --> ("
                  << endKey.first << ", " << endKey.second << ")\n";

        cout << "  [Start Node] x: " << startNode.x
                  << ", y: " << startNode.y
                  << ", psi: " << startNode.psi << "\n";
        cout << "  [End Node]   x: " << endNode.x
                  << ", y: " << endNode.y
                  << ", psi: " << endNode.psi << "\n";

        cout << "  coeffs_x (" << spline.coeffs_x.rows() << "x" << spline.coeffs_x.cols() << "):\n";
        cout << spline.coeffs_x << "\n";

        cout << "  coeffs_y (" << spline.coeffs_y.rows() << "x" << spline.coeffs_y.cols() << "):\n";
        cout << spline.coeffs_y << "\n";

        cout << "----------------------------------------";
    }

    cout << "\n=== ✅ End of splineMap ===\n";
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
    IVector raceline_index_array;
    Vector2d node_pos;

    genNode(nodesPerLayer,
            raceline_index_array,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    // sampling points' info 
    // writeDMapToCSV("inputs/sampling_map.csv", sampling_map);

    Graph edgeList;
    SplineMap splineMap;

    genEdges(nodesPerLayer,
             raceline_index_array,
             edgeList,
             splineMap,
             params.LAT_OFFSET,
             params.LAT_RESOLUTION,
             params.CURVE_THR);

    printSplineMapVerbose(splineMap, nodesPerLayer);
    // edgeList.printGraph();
    
    // visual process 
    visual(edgeList, nodesPerLayer, splineMap);

    return 0;
}