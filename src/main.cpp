#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label);
}

// DMap을 CSV에 작성 
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',') {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();

    // Header
    size_t i = 0;
    for (const auto& [key, _] : map) {
        file << key;
        if (++i != num_cols) file << delimiter;
    }
    file << '\n';

    // Row map
    for (size_t row = 0; row < num_rows; ++row) {
        size_t j = 0;
        for (const auto& [_, col] : map) {
            file << col[row];
            if (++j != num_cols) file << delimiter;
        }
        file << '\n';
    }

    file.close();
}

// Debug용 함수: map의 columns, rows 개수 print  
void map_size(DMap& map) {
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

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

        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;
        int node_idx = 0;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;
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
            node.psi = 0.0;
            node.kappa = 0.0;        
            node.raceline = (node_idx == raceline_index);

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

        // cout << i << "번째 Layer의" << endl;
        // for (size_t i =0; i < node_pos.size(); ++i) {        
        //     cout << i << "번째 Node" << endl;
        //     cout << node_pos[i] << endl;
        // }
                // 각 node의 psi, kappa 계산하는 로직 추가

    }
}
#if 0

VectorXd computeEuclideanDistances(const vector<Vector2d>& path) {
    int N = static_cast<int>(path.size()) - 1;
    VectorXd dists(N);
    for (int i = 0; i < N; ++i) {
        dists(i) = (path[i + 1] - path[i]).norm();
    }
    return dists;
}

unique_ptr<SplineResult> calcSplines(
    const vector<Vector2d>& path,
    const VectorXd *el_lengths = nullptr,
    double psi_s = NAN,
    double psi_e = NAN,
    bool use_dist_scaling = true)
{
    vector<Vector2d> updated_path = path;
    bool closed = (path.front() - path.back()).norm() < 1e-6;

    if (closed) {
        updated_path.push_back(path[1]);  // path[1] 복사해서 끝에 추가
    }

    int no_splines = static_cast<int>(updated_path.size()) - 1;
    VectorXd ds = (el_lengths == nullptr) ? computeEuclideanDistances(updated_path) : *el_lengths;

    MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
    VectorXd b_x = VectorXd::Zero(no_splines * 4);
    VectorXd b_y = VectorXd::Zero(no_splines * 4);

    for (int j = 0; j < no_splines; ++j) {
        int row0 = j * 4;
        M(row0, row0) = 1;
        b_x(row0) = updated_path[j].x();
        b_y(row0) = updated_path[j].y();

        RowVector4d T;
        double d = ds(j);
        T << 1, d, d * d, d * d * d;

        for (int k = 0; k < 4; ++k)
            M(row0 + 1, row0 + k) = T(k);

        b_x(row0 + 1) = updated_path[j + 1].x();
        b_y(row0 + 1) = updated_path[j + 1].y();
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


void genEdges(const vector<vector<Node>> &nodesPerLayer,
              Vector2d &node_pos,
              const double lat_offset,
              const double stepsize_approx,
              const double min_vel = 0.0,
              bool closed = true) {

    if (lat_offset <= 0.0) {
        throw invalid_argument("Too small lateral offset!");
    }

    vector<Vector2d> raceline_cl;
    Vector2d start_point();

    if (closed) {
        raceline_cl.push_back(node_pos.x(), node_pos.y());  // node_pos가 vector<Vector2d>라 가정
    }   
    // node_pos 스플라인 계수 계산
    auto raceline_result = calcSplines(raceline_cl);

    const MatrixXd& x_coeff_r = raceline_result->coeffs_x;
    const MatrixXd& y_coeff_r = raceline_result->coeffs_y;

    for (size_t i = 0; i < nodesPerLayer.size(); ++i) {
        size_t start_layer = i;
        size_t end_layer = i + 1;

        if (end_layer >= nodesPerLayer.size()) {
            if (closed)
                end_layer -= nodesPerLayer.size();
            else
                break;
        }

        for (size_t start_n = 0; start_n < nodesPerLayer[start_layer].size(); ++start_n) {
            const Node& start_node = nodesPerLayer[start_layer][start_n];

            // raceline index 비교를 위해 기준 index 필요
            int ref_start_idx = /* 예: raceline에서의 index 저장 벡터가 필요 */;
            int ref_end_idx = /* 마찬가지로 end_layer에서의 raceline index */;

            int end_n_ref = ref_end_idx + start_n - ref_start_idx;

            // 거리 계산
            Vector2d d_start(start_node.x, start_node.y);
            const auto& end_layer_nodes = nodesPerLayer[end_layer];
            int clipped_idx = clamp(end_n_ref, 0, (int)end_layer_nodes.size() - 1);
            Vector2d d_end(end_layer_nodes[clipped_idx].x, end_layer_nodes[clipped_idx].y);

            double dist = (d_end - d_start).norm();
            int lat_steps = static_cast<int>(round(dist * lat_offset / start_node.lat_resolution));

            for (int end_n = max(0, end_n_ref - lat_steps); end_n <= min((int)end_layer_nodes.size() - 1, end_n_ref + lat_steps); ++end_n) {
                const Node& end_node = end_layer_nodes[end_n];

                MatrixXd x_coeff, y_coeff;

                if (start_node.raceline && end_node.raceline) {
                    x_coeff = x_coeff_r.row(start_layer);
                    y_coeff = y_coeff_r.row(start_layer);
                } else {
                    vector<Vector2d> path = { Vector2d(start_node.x, start_node.y), Vector2d(end_node.x, end_node.y) };
                    auto result = calcSplines(path, nullptr, start_node.psi, end_node.psi);
                    x_coeff = result->coeffs_x;
                    y_coeff = result->coeffs_y;
                }

                // addEdge(start_layer, start_n, end_layer, end_n, x_coeff, y_coeff);
            }
        }
    }
}

#endif




int main() {
    #if 1
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
    #endif

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