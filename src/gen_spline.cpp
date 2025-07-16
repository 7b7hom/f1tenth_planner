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

    size_t num_cols = map.size();   // 열 개수 = map의 key 개수
    size_t num_rows = map.begin()->second.size();   // 행 개수 = 첫 번째 key의 값의 길이

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
    // cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

// Dvector를 Map 구조로 추가(연산)
void addDVectorToMap(DMap &map,
                     string attr,
                     const IVector *idx_array = nullptr) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();  // 지정X -> 전체 데이터 대상
    } 
    else {
        len = idx_array->size();    // 지정한 인덱스 수만큼
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

// samplePoint : 레이어가 생성될 위치 정하는 것. (결과적으로 레이어 간격 의미)
void samplePointsFromRaceline(const DVector& kappa,     // 곡률
                              const DVector& dist,      // 점 사이 거리
                              double d_curve,           // 곡선 구간 샘플링 간격
                              double d_straight,        // 직선 구간 샘플링 간격
                              double curve_th,          // 곡선 판단 기준 곡률
                              IVector& idx_array) {     // 결과 (샘플링 인덱스) 저장할 벡터

    const size_t n = kappa.size();  // 전체 경로 길이 (전체 경로에서 몇 개의 점을 갖는지)
    double cur_dist = 0.0;          // 지금까지 이동한 거리
    double next_dist = 0.0;         // 다음 샘플링 지점 결정하는 기준 거리
    double next_dist_min = 0.0;     // 곡선 구간일 경우 : 다음 샘플링까지 최소 간격 확보용

    for (size_t i = 0; i < n; ++i) {    // 경로의 각 점에 대해 반복

        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {    // fabs() : float absolute value (절댓값 반환)
            next_dist = cur_dist;
        }

        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) { // 현재 위치가 next_dist 넘었다면 샘플링 지점 도달!
            idx_array.push_back(static_cast<int>(i));   // 명시적 형변환

            if (fabs(kappa[i]) < curve_th) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve; // 너무 가까운 곳은 샘플 중복해서 찍지 않도록 해줌
        }

        cur_dist += dist[i];
    }

    //  for (size_t i=0; i < idx_array.size(); ++i) 
    //      cout << idx_array[i] << endl;
    //  cout << "size: " << idx_array.size() << endl;
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
    psi[i] = normalizeAngle(atan2(dy, dx) - M_PI_2);

    // cout << i<< ": " << psi[i] << endl;

    }
    // cout << psi.size() << endl;

}

void genNode(NodeMap& nodesPerLayer,        // 각 레이어에 생성된 노드 저장하는 2차원 벡터
            IVector& raceline_index_array,  // 각 레이어에서 레이싱라인이 위치한 노드의 인덱스
            const double veh_width,         // 차량의 너비
            float lat_resolution) {         // 노드 간 lateral 간격 (옆 방향 노드들 간격, m 단위)
    
    const size_t N = sampling_map[__alpha].size();  // 총 레이어 수
    Vector2d node_pos;
    nodesPerLayer.resize(N);    // N개 레이어 기준, nodesPerLayer 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        node.layer_idx = i; 
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index); // 레이스라인 노드 저장
        
        // cout << "layer 길이" << (sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width/2)<< endl;
        // cout << "layer 내에서 raceline index:" << raceline_index << endl;


        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);    // 레이어의 기준이 되는 점 (기준선 위 한 점)
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);  // ref point에서 옆 방향을 알려주는 벡터
        
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지 (lat_resolution씩 증가시키며 노드 생성)
        int node_idx = 0;
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;  // num_nodes : 좌우 총 가능한 노드 수
        nodesPerLayer[i].resize(num_nodes); 

        // cout << i << "번째 layer의 node 개수는 " << num_nodes << endl;

        // node별 loop (노드 생성 시작)
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

            // psi 보간 처리
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
}       // 저장 결과 : nodesPerLayer[i][j] (i번째 레이어에서 j번째 lateral 위치의 노드)
        //           raceline_index_array[i] (i번째 레이어에서 레이싱라인이 위치한 인덱스)
        //           노드 정보

void plotHeading(const DVector &x, const DVector &y, const DVector &psi, double scale = 0.5) {
    double dx, dy;
    double theta, arrow_len;
    double angle;
    double x_arrow1, y_arrow1;
    double x_arrow2, y_arrow2;

    for (size_t i = 0; i < x.size(); ++i) {
        dx = scale * cos(psi[i] + M_PI_2);
        dy = scale * sin(psi[i] + M_PI_2);
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        plt::plot(x_line, y_line, {{"color", "green"}});

        #if 1 // 화살촉 그리기
        theta = atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);
        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "green"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "green"}});
        #endif
    }
}

// NodeMap에 저장된 모든 노드들을 보라색 점으로 플로팅하고, 각 노드의 헤딩을 화살표로 시각화
void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            node_x.push_back(node.x);
            node_y.push_back(node.y);
        }
    }
    plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}, {"label", "Nodes"}});
}



// -------------------------------------------
// gen_spline --------------------------------
// -------------------------------------------

// 경로 행렬의 각 점 사이의 유클리드 거리 계산 함수
VectorXd computeEuclideanDistances(const MatrixXd& path) {
    int N = path.rows() - 1;    // 구간 개수는 점 개수 - 1
    VectorXd dists(N);          // 거리 저장 벡터
    for (int i = 0; i < N; ++i) {
        dists(i) = (path.row(i + 1) - path.row(i)).norm();  // 두 점 사이 거리 계산 (유클리드 노름)
    }
    return dists;
}

// only startNode, endNode
SplineResult calcSplines(const Node& startNode, const Node& endNode) {
    // fixed M matrix assuming d = 1
    MatrixXd M(4, 4);
    VectorXd b_x(4), b_y(4);

    // heading 방향 90도 회전
    double psi_s = startNode.psi + M_PI_2;
    double psi_e = endNode.psi + M_PI_2;

    // 거리 정규화
    // d = 1 fixed
    M << 1, 0, 0, 0,    // x(0) = a_0
         1, 1, 1, 1,    // x(1) = a_0 + a_1 + a_2 + a_3
         0, 1, 0, 0,    // x'(0) = a_1
         0, 1, 2, 3;    // x'(1) = a_1 + 2a_2 + 3a_3

    b_x << startNode.x,
           endNode.x,
           cos(psi_s),
           cos(psi_e);

    b_y << startNode.y,
           endNode.y,
           sin(psi_s),
           sin(psi_e);

    VectorXd coeffs_x = M.colPivHouseholderQr().solve(b_x);
    VectorXd coeffs_y = M.colPivHouseholderQr().solve(b_y);

    // transpose
    MatrixXd coeffs_x_trans = coeffs_x.transpose();  // (1x4)
    MatrixXd coeffs_y_trans = coeffs_y.transpose();  // (1x4)

    SplineResult result;
    result.coeffs_x = coeffs_x_trans;
    result.coeffs_y = coeffs_y_trans;

    return result;
}


// ------------------ genEdge ------------------

Vector2d computeSplinePosition(const RowVector4d& coeff_x, const RowVector4d& coeff_y, double t) {
    double t2 = t * t;
    double t3 = t2 * t;
    double x = coeff_x(0) + coeff_x(1) * t + coeff_x(2) * t2 + coeff_x(3) * t3;
    double y = coeff_y(0) + coeff_y(1) * t + coeff_y(2) * t2 + coeff_y(3) * t3;
    return Vector2d(x, y);
}


void genEdge(Graph& graph, 
    const NodeMap& nodesPerLayer, 
    const Offline_Params& params,
    const IVector& raceline_index_array,
    bool closed = true) {


        // closed -> 마지막과 0번 레이어 연결, open -> break
        for (size_t layer = 0; layer < nodesPerLayer.size(); ++layer) {
            
            const size_t start_layer = layer;
            const size_t end_layer = (layer + 1) % nodesPerLayer.size();
        
            if (params.CLOSURE_DETECTION_DIST < 1e-6 && end_layer == 0) break;

            // 1. 각 레이어의 레이스라인 기준 인덱스
            const int start_race_idx = raceline_index_array[start_layer];
            const int end_race_idx   = raceline_index_array[end_layer];

            // 2. 시작 레이어의 각 노드에 대해 반복
            const auto& start_layer_nodes = nodesPerLayer[start_layer];
            const auto& end_layer_nodes   = nodesPerLayer[end_layer];

            for (size_t start_idx = 0; start_idx < start_layer_nodes.size(); ++start_idx) {
                const Node& startNode = start_layer_nodes[start_idx];

            // 3. 레이스라인에서 같은 위치에 있는 end 레이어의 참조 노드 인덱스 계산
            int offset_from_raceline = static_cast<int>(start_idx) - start_race_idx;
            int aligned_end_idx = end_race_idx + offset_from_raceline;

            // 4. 범위 클램프 (value, min, max) : value가 [min, max] 사이에 있도록 제한 (min보다 작으면 min, max보다 크면 max 반환)
            aligned_end_idx = std::clamp(aligned_end_idx, 0, static_cast<int>(end_layer_nodes.size() - 1));

            // 5. 거리 계산을 위한 노드 좌표 가져오기
            const Node& ref_end_node = end_layer_nodes[aligned_end_idx];

            // 6. 거리 계산을 위한 좌표 행렬 생성
            MatrixXd spline_path(2, 2);
            spline_path << startNode.x, startNode.y,
                        ref_end_node.x, ref_end_node.y;


            VectorXd el_lengths = computeEuclideanDistances(spline_path);
            double dist = el_lengths(0);

            // 커브면 더 많이 연결
            double factor = (startNode.kappa > params.CURVE_THR) ? 2.0 : 1.0;
            int lat_steps = round(factor * dist * params.LAT_OFFSET / params.LAT_RESOLUTION);

            for (int destIdx = std::max(0, aligned_end_idx - lat_steps);
                destIdx <= std::min(static_cast<int>(end_layer_nodes.size() - 1), aligned_end_idx + lat_steps);
                ++destIdx) {

                const Node& endNode = end_layer_nodes[destIdx];

                // 스플라인 계산
                auto result = calcSplines(startNode, endNode);
                const MatrixXd& x_coeffs = result.coeffs_x;
                const MatrixXd& y_coeffs = result.coeffs_y;

                // DEBUG

                // double d = std::sqrt(std::pow(endNode.x - startNode.x, 2) +
                //      std::pow(endNode.y - startNode.y, 2));

                // for (double t = 0; t <= d; t += 0.2) {
                //     Vector2d pt = computeSplinePosition(x_coeffs.row(0), y_coeffs.row(0), t);
                //     std::cout << "  - t=" << t << " → (" << pt.x() << ", " << pt.y() << ")\n";
                // }
                // std::cout << "start: (" << startNode.x << ", " << startNode.y << "), "
                //         << "end: (" << endNode.x << ", " << endNode.y << ")" << std::endl;

                // std::cout << "coeff_x: " << x_coeffs << std::endl;
                // std::cout << "coeff_y: " << y_coeffs << std::endl;

                // for (int k = 0; k <= 10; ++k) {
                //     double t = static_cast<double>(k) / 10.0;
                //     double x = x_coeffs(0) + x_coeffs(1) * t + x_coeffs(2) * t * t + x_coeffs(3) * t * t * t;
                //     double y = y_coeffs(0) + y_cofeffs(1) * t + y_coeffs(2) * t * t + y_coeffs(3) * t * t * t;
                //     std::cout << "  - point(" << t << "): (" << x << ", " << y << ")\n";
                // }
                

                // 그래프에 엣지 추가
                ITuple src_key(start_layer, startNode.node_idx);
                // std::cout << "  → Adding edge: (" << std::get<0>(src_key) << "," << std::get<1>(src_key)
                //         << ") → " << endNode.node_idx << "\n";
                graph.addEdge(src_key, endNode.node_idx);

            }
        }  
    }
}

void visual(const NodeMap& nodesPerLayer, Graph& graph, const Offline_Params& params) {
    plt::clf();

    // 트랙 경계선
    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // 레이싱 라인 및 샘플링된 포인트
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}, {"label", "Raceline"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}, {"label", "Sampled Raceline"}});
    plotHeading(sampling_map[__x_raceline], sampling_map[__y_raceline], sampling_map[__psi]);

    plotHeading(nodesPerLayer);
    
    DVector spline_x_pts; 
    DVector spline_y_pts;

    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& current_node : layer_nodes) {
            ITuple src_key(current_node.layer_idx, current_node.node_idx);
            IVector child_nodes_idx;

            try {
                graph.getChildIdx(src_key, child_nodes_idx);
            } catch (const std::runtime_error& e) {
                continue;
            }
            
            for (int dest_node_idx : child_nodes_idx) {

                spline_x_pts.clear(); 
                spline_y_pts.clear(); 
                
                size_t next_layer_idx = (current_node.layer_idx + 1) % nodesPerLayer.size();

                if (dest_node_idx < 0 || dest_node_idx >= nodesPerLayer[next_layer_idx].size()) {
                    // cout << "Warning: Invalid dest_node_idx " << dest_node_idx << " for layer " << next_layer_idx << endl; // 이 라인도 extended character 오류의 원인이 될 수 있습니다.
                    continue;
                }
                const Node& next_node = nodesPerLayer[next_layer_idx][dest_node_idx];

                MatrixXd spline_path(2, 2);
                spline_path << current_node.x, current_node.y,
                               next_node.x, next_node.y;
                
                double psi_s = current_node.psi;
                double psi_e = next_node.psi;

                VectorXd el_lengths(1);
                el_lengths(0) = (spline_path.row(1) - spline_path.row(0)).norm();

                SplineResult res;
                try {
                    res = calcSplines(current_node, next_node);
                } catch (const std::exception& e) {
                    continue;
                }

                const int num_spline_segments = 10; 
                for (int k = 0; k <= num_spline_segments; ++k) {
                    double t_eval = static_cast<double>(k) / num_spline_segments;
                    
                    Vector2d pos = computeSplinePosition(res.coeffs_x.row(0), res.coeffs_y.row(0), t_eval);
                    
                    spline_x_pts.push_back(pos.x());
                    spline_y_pts.push_back(pos.y());
                }

                plt::plot(spline_x_pts, spline_y_pts, {{"color", "green"}, {"linewidth", "1"}}); // {"label", "Valid Splines"}
            }
        }
    }

    plt::title("Track and Planned Graph");
    plt::grid(true);
    plt::axis("equal");
    plt::legend();
    plt::show();
}



int main() {
    IVector idx_sampling;
    Offline_Params params;

    string map_file_in = "../inputs/gtpl_levine.csv";
    string map_file_out = "../inputs/gtpl_levine_out.csv";

    // global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    readDMapFromCSV(map_file_in, gtpl_map);
    std::cout << "CSV loaded, columns: " << gtpl_map.size() << std::endl;

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
    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);

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

    // 노드 그리드 생성
    NodeMap nodesPerLayer;
    IVector raceline_index_array;
    Vector2d node_pos;

    genNode(nodesPerLayer,
            raceline_index_array,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    //Graph edgeList;

    MatrixXd path_xy(sampling_map[__x_raceline].size(), 2);
    for (size_t i = 0; i < sampling_map[__x_raceline].size(); ++i) {
        path_xy(i, 0) = sampling_map[__x_raceline][i];
        path_xy(i, 1) = sampling_map[__y_raceline][i];
    }

    double psi_s = sampling_map[__psi][0];
    double psi_e = sampling_map[__psi].back();

    // 최종 그래프 생성
    Graph myGraph;
    genEdge(myGraph,
            nodesPerLayer,
            params,
            raceline_index_array);
    
    // 시각화
    visual(nodesPerLayer, myGraph, params);

    return 0;
}