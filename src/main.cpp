#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames(); // CSV 파일의 모든 열 이름(예: "x_ref", "y_ref", "kappa" 등)을 가져와 labels 벡터에 저장

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label); // 해당 열 이름에 해당하는 모든 데이터를 double 타입의 벡터로 가져오고 map에 label을 키로하여 저장
}

// DMap을 CSV에 작성 
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',') {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    // 열의 개수(맵의 키 개수)와 행의 개수(각 벡터의 크기)를 가져옴.
    size_t num_cols = map.size(); // key 부분 접근
    size_t num_rows = map.begin()->second.size(); // value 부분(->second) 접근

    // Header 맵의 모든 키(열 이름)를 지정된 delimiter(구분자)로 구분하여 파일에 씀.
    // -- 헤더(열 이름) 작성 -- 
    size_t i = 0;
    // map의 모든 키-값 쌍에 대해 반복
    for (const auto& [key, _] : map) {
        file << key;
        if (++i != num_cols) file << delimiter; // 마지막 열이 아니면 구분자 씀.
    }
    file << '\n';

    // Row map
    // -- 데이터(행) 작성 -- 
    // 각 행에 대해 반복
    for (size_t row = 0; row < num_rows; ++row) {
        size_t j = 0;
        for (const auto& [_, col] : map) {
            file << col[row]; // 현재 열(col)의 현재 행(row)에 해당하는 데이터를 파일에 씀.
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
// map에 이미 있는 데이터를 기반으로 새로운 기하학적 속성(예: 트랙 경계선, 레이싱 라인)의 좌표를 계산하여 map에 추가
void addDVectorToMap(DMap &map,
                     string attr, // 계산할 속성의 이름(예: bound_r, bound_l, raceline, delta_s 등)
                     const IVector *idx_array = nullptr) {
    size_t len; // 새로 생성될 벡터의 길이(=> DVector (즉, std::vector<double>)에 몇 개의 double 값들이 들어갈 것인가)
    if (idx_array == nullptr) { // 이 함수가 gtpl_map (원본 트랙 데이터)에 bound_r, bound_l, raceline, delta_s와 같은 새 데이터를 추가할 때 사용
        len = map[__x_ref].size(); // 원본 트랙의 전체 점 개수(행 개수)를 len으로 사용
    } 
    else { // 이 함수가 sampling_map (샘플링된 트랙 데이터)에 delta_s 등을 추가할 때 사용 
        len = idx_array->size(); // 샘플링된 포인트의 개수를 len으로 사용
    }
    // cout << "attr: "<< attr << " / len:" << len << endl;

    DVector x_out(len), y_out(len); // 계산된 X, Y 좌표를 저장할 벡터 (len 크기로 초기화)
    string x_label = "x_" + attr; // X좌표 레이블 (예: "x_bound_r")
    string y_label = "y_" + attr; // Y좌표 레이블 (예: "y_bound_r")
    
    // attr이 "bound_r" (오른쪽 경계선)인 경우
    if (!attr.compare("bound_r")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            // 오른쪽 경계의 X좌표 = 기준선 X좌표 + 기준선 법선 벡터 X성분 * 오른쪽 폭
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__width_right][i];
            // 오른쪽 경계의 Y좌표 = 기준선 Y좌표 + 기준선 법선 벡터 Y성분 * 오른쪽 폭
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__width_right][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out; // 계산된 X좌표 벡터를 맵에 "x_bound_r" 키로 저장
        map[y_label] = y_out; // 계산된 Y좌표 벡터를 맵에 "y_bound_r" 키로 저장
    }
    // attr이 "bound_l" (왼쪽 경계선)인 경우
    else if (!attr.compare("bound_l")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            // 왼쪽 경계선의 X좌표 = 기준선 X좌표 - 기준선 법선 벡터 X성분 * 왼쪽 폭
            x_out[i] = map[__x_ref][i] - map[__x_normvec][i] * map[__width_left][i];
            // 왼쪽 경계선의 Y좌표 = 기준선 Y좌표 - 기준선 법선 벡터 Y성분 * 왼쪽 폭
            y_out[i] = map[__y_ref][i] - map[__y_normvec][i] * map[__width_left][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out; // 계산된 X좌표 벡터를 맵에 "x_bound_l" 키로 저장
        map[y_label] = y_out; // 계산된 Y좌표 벡터를 맵에 "y_bound_l" 키로 저장
    }
    // attr이 "raceline"인 경우
    else if (!attr.compare("raceline")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            // 레이싱 라인 X좌표 = 기준선 X좌표 + 기준선 법선 벡터 X성분 * 알파 (알파는 레이싱 라인의 횡방향 오프셋)
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
    // attr이 "delta_s" (연속된 점 사이의 호 길이 변화량)인 경우
    else if (!attr.compare("delta_s")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len - 1; ++i) {
            // delta_s[i] = 다음 레이싱 궤적의 s 값 - 현재 레이싱 궤적의 s 값
            x_out[i] = map[__s_racetraj][i+1] - map[__s_racetraj][i]; // 마지막 원소는 0
            // (기본적으로 Dvector x_out(len) 선언 시 len개의 원소를 기본값 0.0으로 초기화하기 때문에 
            // len-2번 돌아서 마지막 값이 할당되지 않아도 마지막 원소는 0.0으로 유지)
        }
        map[attr] = x_out; // 계산된 delta_s 벡터를 맵에 "delta_s" 키로 저장
    }

    // map_size(map);
}

// 레이싱 라인에서 경로 계획을 위한 layer가 될 지점들을 샘플링
// 곡선 구간에서는 더 촘촘하게, 직선 구간에서는 느슨하게 샘플링하여 효율성을 높임.
// 이미 존재하는 촘촘한 원본 데이터에서, 트랙의 특성(직선/곡선)에 따라 필요한 만큼의 밀도로 포인트들을 '선별적'으로 추출해냄.
void samplePointsFromRaceline(const DVector& kappa, // 레이싱 라인의 각 지점에서의 curvature 벡터
                              const DVector& dist, // 점과 점 사이 길이 벡터
                              double d_curve, // 곡선 구간에서 샘플링할 최소 거리 간격
                              double d_straight, // 직선 구간에서 샘플링할 최소 거리 간격
                              double curve_th, //  곡선으로 간주할 곡률 임계값, 이 값보다 |kappa|가 크면 곡선으로 간주
                              IVector& idx_array) { // 샘플링된 지점들의 원본 레이싱 라인에서의 인덱스가 저장될 벡터

    const size_t n = kappa.size(); // 레이싱 라인 전체 포인트의 개수
    double cur_dist = 0.0; // 현재까지 레이싱 라인을 따라 누적된 거리
    double next_dist = 0.0; // 다음 샘플링 지점까지의 목표 거리
    double next_dist_min = 0.0; // 곡선 구간에서 최소 간격을 보장하기 위한 거리

    for (size_t i = 0; i < n; ++i) {

        // 곡선이면 최소 거리 갱신
        // cur_dist + dist[i] : 현재 포인트까지의 누적 거리
        // -> 현재 위치(i번째 포인트)가 다음 샘플링을 고려할 수 있는 최소 거리(next_dist_min)를 넘어섰는지 확인
        // fabs(kappa[i]) > curve_th -> 곡선 구간에 진입했는지 확인
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist; // 다음 샘플링 지점 목표를 현재 위치로 재설정 (곡선 시작 지점에서 바로 샘플링하기 위함)
        }

        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) {
            idx_array.push_back(static_cast<int>(i)); // 현재 포인트의 인덱스를 샘플링 배열에 추가

            if (fabs(kappa[i]) < curve_th) {  // 직선 구간
                next_dist += d_straight; // 다음 샘플링 목표 거리를 직선 간격만큼 늘림.
            } else {  // 곡선 구간
                next_dist += d_curve; // 다음 샘플링 목표 거리를 곡선 간격만큼 늘림.
            }

            // next_dist_min은 곡선 구간에서 최소 d_curve 간격을 유지하기 위한 값(곡선 구간에서 너무 드문드문 샘플링되지 않게 하기 위함.)
            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist[i]; // 현재까지 누적된 거리 업데이트
    }

    // for (size_t i=0; i < idx_array.size(); ++i) 
    //     cout << idx_array[i] << endl;
    // cout << "size: " << idx_array.size() << endl;
}

// 주어진 각도를 -PI에서 PI 사이의 값으로 정규화 (예: 360도는 0도, 270도는 -90도, 450도는 90도)
double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

// 주어진 X, Y 좌표 벡터를 기반으로 각 지점에서의 헤딩(Heading, 진행 방향 각도)을 계산
void calcHeading(DVector &x_raceline, // 라인의 X좌표 벡터
                 DVector &y_raceline, // 라인의 Y좌표 벡터
                 DVector &psi) { // 계산된 헤딩 각도(라디안)가 저장될 벡터

    size_t N = x_raceline.size(); // 포인트 개수
    psi.resize(N); // psi 벡터의 크기를 포인트 개수와 같게 조정

    // 닫힌 회로 가정. 예외 처리 필요
    double dx, dy;
    for (size_t i = 0; i < N; ++i) {
        
        if (i != N -1) { // 마지막 포인트가 아닐 시
            dx = x_raceline[i+1] - x_raceline[i]; // 현재 포인트와 다음 포인트의 X 변화량
            dy = y_raceline[i+1] - y_raceline[i]; // 현재 포인트와 다음 포인트의 Y 변화량
        } else { // 마지막 포인트라면 (닫힌 회로, 즉 트랙이 한 바퀴 도는 형태를 가정)
            dx = x_raceline[0] - x_raceline[N - 1]; // 마지막 포인트와 첫 번째 포인트의 X 변화량
            dy = y_raceline[0] - y_raceline[N - 1]; // 마지막 포인트와 첫 번째 포인트의 Y 변화량
        }
    // atan2(dy, dx): (dy, dx) 벡터가 x축과 이루는 각도 반환
    // - M_PI_2: 각도를 트랙의 진행 방향(Y축 방향)에 맞게 조정
    psi[i] = atan2(dy, dx) - M_PI_2;
        
    normalizeAngle(psi[i]); // 계산 값 정규화

    }
    // cout << i<< ": " << psi[i] << endl;
    // cout << psi.size() << endl;

}

// 샘플링된 레이어마다 경로 계획을 위한 Node(차량이 횡방향으로 이동 가능한 위치들) 생성
void genNode(NodeMap& nodesPerLayer, // 생성된 node들이 저장될 NodeMap(vector<vector<Node>>)
            const double veh_width, // 차량의 폭
            float lat_resolution) { // 횡방향(트랙 폭 방향) node 간의 간격

    const size_t N = sampling_map[__alpha].size(); // 샘플링된 레이어의 개수
    IVector raceline_index_array; // 각 레이어에서 레이싱 라인 node의 횡방향 인덱스를 저장할 벡터
    Vector2d node_pos; // node의 (x, y) 좌표를 임시로 저장할 변수
    nodesPerLayer.resize(N); // NodeMap의 크기를 레이어 개수 N에 맞게 조정

    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        node.layer_idx = i;

        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        // sampling_map[__width_left][i]: i번째 레이어에서 레퍼런스 라인(기준선)으로부터 트랙의 왼쪽 경계까지의 거리
        //                                -> 트랙의 가장 왼쪽 경계를 0으로 하는 새로운 횡방향 좌표계로 변환
        // + sampling_map[__alpha][i]: 레퍼런스 라인으로부터 레이싱 라인(차량의 중심이 주행하는 경로)까지의 횡방향 오프셋
        //                             -> 트랙의 가장 왼쪽 경계에서 레이싱 라인까지의 횡방향 오프셋으로 바뀜.
        // - veh_width / 2 -> 차량의 물리적 크기를 고려, 차량의 중심이 물리적으로 위치할 수 있는 가장 왼쪽 한계점을 0으로 하는 새로운 좌표계로 기준점 이동
        // 최종 raceline_index : 차량 중심이 물리적으로 위치할 수 있는 가장 왼쪽 한계점"을 0번 인덱스로 하는 그리드 상에서, 
        //                       레이싱 라인(차량 중심)이 몇 번째 인덱스에 해당하는지
        int raceline_index = floor((sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width / 2) / lat_resolution);
        raceline_index_array.push_back(raceline_index);
        
        // cout << "layer 길이" << (sampling_map[__width_left][i] + sampling_map[__alpha][i] - veh_width/2)<< endl;
        // cout << "layer 내에서 raceline index:" << raceline_index << endl;
        // cout << "-----" << endl;

        // 현재 레이어의 기준점 (ref_xy)과 법선 벡터 (norm_vec)를 가져옴.
        Vector2d ref_xy(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);
        
        // start_alpha: 현재 레이어에서 노드를 생성할 횡방향 오프셋의 시작 값
        double start_alpha = sampling_map[__alpha][i] - raceline_index * lat_resolution;
        int node_idx = 0; // 현재 레이어 내에서 노드의 인덱스
        // num_nodes: 현재 레이어 내에서 생성될 총 노드 개수
        // (전체 트랙 폭 - 차량 폭) / node 간격 + 1
        int num_nodes = (sampling_map[__width_right][i] + sampling_map[__width_left][i] - veh_width) / lat_resolution + 1;
        nodesPerLayer[i].resize(num_nodes); // 현재 레이어에 해당하는 벡터의 크기를 노드 개수에 맞게 조정(미리 필요한 공간 확보(동적x))

        // cout << i << "번째 layer의 node 개수는 " << num_nodes << endl;
        // -- 현재 레이어 내에서 각 노드를 생성하는 루프 --
        // node별 loop 
        // 오른쪽 경계에서 차량 폭의 절반을 뺀 값까지 노드를 생성
        for (double alpha = start_alpha; alpha <= sampling_map[__width_right][i] - veh_width / 2 ; alpha+=lat_resolution) {
            // node_alphas.push_back(alpha);
            // node의 좌표 계산: 기준점 + 횡방향 오프셋 * 법선 벡터
            node_pos = ref_xy + alpha * norm_vec;

            // node의 layer내의 인덱스 계산.
            node.node_idx = node_idx;
            node.x = node_pos.x();
            node.y = node_pos.y();
            node.psi = 0.0;
            node.kappa = 0.0;        
            node.raceline = (node_idx == raceline_index);

            // -- 노드의 헤딩(psi) 계산 (보간) --
            double psi_interp; // 보간된 psi 값
            if (node_idx < raceline_index) { // 레이싱 라인 노드보다 왼쪽에 있는 노드들
                // 왼쪽 경계선 psi와 레이싱 라인 psi 사이 보간
                if (abs(sampling_map[__psi_bound_l][i] - sampling_map[__psi][i]) >= M_PI) 
                {   
                    double bl = sampling_map[__psi_bound_l][i] + 2 * M_PI * (sampling_map[__psi_bound_l][i] < 0);
                    double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
                    psi_interp = bl + (p - bl) * node_idx / raceline_index;
                                       
                }
                else {
                    psi_interp = sampling_map[__psi_bound_l][i] + (sampling_map[__psi][i] - sampling_map[__psi_bound_l][i]) * (node_idx+1) / raceline_index;
                }
                node.psi = normalizeAngle(psi_interp); // 보간된 psi 정규화
            }
            else if (node_idx == raceline_index) { // 레이싱 라인 노드
                psi_interp = sampling_map[__psi][i]; // 레이싱 라인의 psi 값을 그대로 사용
                node.psi = psi_interp;
            }
            else { // 레이싱 라인 노드보다 오른쪽에 있는 노드들
                // 오른쪽 경계선 psi와 레이싱 라인 psi 사이 보간
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

// 주어진 (x, y) 좌표에서 해당 psi (헤딩) 방향을 화살표로 시각화
void plotHeading(const DVector &x, // 플로팅할 지점들의 X, Y 좌표 벡터
                 const DVector &y,
                 const DVector &psi, // 각 지점의 헤딩 각도 벡터
                 double scale = 0.5) // 화살표의 길이 스케일
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
        plt::plot(x_line, y_line, {{"color", "green"}});

        #if 1
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
        #if 0
        for (size_t i = 0; i < x.size(); ++i) {
        ostringstream label;
        label.precision(2);
        label << fixed << "(" << x[i] << ", " << y[i] << ")\nψ=" << psi[i];

        plt::text(x[i], y[i], label.str());
    }
        #endif
}

// NodeMap에 저장된 모든 노드들을 보라색 점으로 플로팅하고, 각 노드의 헤딩을 화살표로 시각화
void plotHeading(const NodeMap& nodesPerLayer, double scale = 0.5) {
    DVector x_line, y_line;
    DVector node_x, node_y;
    for (const auto& layer_nodes : nodesPerLayer) {
        for (const auto& node : layer_nodes) {
            double dx = scale * cos(node.psi + M_PI_2);
            double dy = scale * sin(node.psi + M_PI_2);

            node_x.push_back(node.x);
            node_y.push_back(node.y);
            plt::scatter(node_x, node_y, 15.0, {{"color", "purple"}});

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
        
    }
}

// 트랙의 경계, 레이싱 라인, 샘플링된 포인트, 그리고 생성된 노드들을 한 화면에 종합적으로 시각화
void visual(const NodeMap& nodesPerLayer) {
    plt::clf();

    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // plt::plot(gtpl_map[__x_ref], gtpl_map[__y_ref], {{"color", "blue"}});
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}});

    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}});

    plotHeading(sampling_map[__x_raceline],
                sampling_map[__y_raceline],
                sampling_map[__psi]);

    // plotHeading(sampling_map[__x_bound_l],
    //             sampling_map[__y_bound_l],
    //             psi_bound_l);

    // plotHeading(sampling_map[__x_bound_r],
    //             sampling_map[__y_bound_r],
    //             psi_bound_r);


    // 노드마다 psi확인할 수 있는 용도 
    plotHeading(nodesPerLayer);

    plt::title("Track");
    plt::grid(true);
	plt::axis("equal");
	plt::show();  
}


int main() {
    #if 0
    IVector idx_sampling; // 샘플링된 인덱스를 저장할 벡터
    Offline_Params params; // 경로 계획에 필요한 설정 값들을 담는 구조체

    string map_file_in = "inputs/gtpl_levine.csv"; // 입력 CSV 파일 경로
    string map_file_out = "inputs/gtpl_levine_out.csv"; // 출력 CSV 파일 경로

    // 1. global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    readDMapFromCSV(map_file_in, gtpl_map);

    // 2. gtpl_map에 오른쪽 경계, 왼쪽 경계, 레이싱 라인, delta_s 데이터 계산하여 추가
    addDVectorToMap(gtpl_map, "bound_r");
    addDVectorToMap(gtpl_map, "bound_l");
    addDVectorToMap(gtpl_map, "raceline");
    addDVectorToMap(gtpl_map, "delta_s");

    writeDMapToCSV(map_file_out, gtpl_map); // 확장된 gtpl_map을 CSV로 저장

    // 3. layer 간격을 위한 raceline points sampling 
    samplePointsFromRaceline(gtpl_map[__kappa],
                             gtpl_map[__delta_s],
                             params.LON_CURVE_STEP,
                             params.LON_STRAIGHT_STEP,
                             params.CURVE_THR,
                             idx_sampling);

    // cout << "idx size:" << idx_sampling.size() << endl;
    
    // 4. 샘플링된 인덱스를 사용하여 gtpl_map에서 데이터를 복사하여 sampling_map 채우기
    for (const auto& [key, vec] : gtpl_map) {
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < vec.size()) { // 유효한 인덱스인지 확인(실제 데이터에 적용하기 전 안전 검검)
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }
    // writeDMapToCSV("inputs/sampling_map", sampling_map);
    // map_size(sampling_map); // (51, 3)

    // 5. sampling_map에 delta_s 다시 계산 (샘플링된 데이터에 맞춰)
    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);
    
    // map_size(sampling_map); // (51, 4)

    // 추후 저장될 예정 
    // 6. sampling_map의 레이싱 라인, 왼쪽/오른쪽 경계의 헤딩(psi) 계산
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
    
    NodeMap nodesPerLayer; // 생성된 노드들을 저장할 NodeMap

    // 7. 노드 생성: 각 레이어에 차량 폭과 lat_resolution 고려하여 노드 그리드 생성
    genNode(nodesPerLayer,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    // sampling points' info 
    // writeDMapToCSV("inputs/sampling_map.csv", sampling_map);
    
    // 8. visual process 
    visual(nodesPerLayer);
    #endif

    // Graph sample code 
    Graph directedGraph;
    ITuple t1(0, 0);
    ITuple t2(0, 1);
    ITuple t3(1, 0);
    ITuple t4(0, 2);
    // spline 생성 후 edge로 집어넣음.
    directedGraph.addEdge(t1, 3);// push_back이라서 sorting은 되지 않음. 
    directedGraph.addEdge(t1, 4);
    directedGraph.addEdge(t1, 5);
    directedGraph.addEdge(t3, 2);
    directedGraph.addEdge(t4, 3);
    // 실제 로직은 node idx가 작은 순서대로 그래프가 그러질 예정이라 괜찮을 듯.
    // grpah 전체 print 
    cout << "---처음 Graph---" << endl;
    directedGraph.printGraph();

    IVector child1;
    // t1 노드의 뒤로 연결된(child) node들을 뽑아온다.
    // directedGraph.getChildIdx(t1, child1);

    // Error Index(child가 없는 경우 runtime_error)
    // directedGraph.getChildIdx(t2, child);
    
    // (0, n)이라는 임의의 노드 n이 (1, 3)을 들고 있는 경우 해당 list에서 3번 삭제
    // 이후 child 노드가 그 뒤로 연결된 spline이 없는 경우 parent의 adjList에서 child 노드를 삭제하기 위하여 필요함.
    vector<ITuple> parent; 
    // 0번째 layer에서 1번째 layer의 3번 노드를 들고 있는 노드가 있는지 찾고 parent로 반환함.
    directedGraph.getParentNode(0, 3, parent); 
    
    cout << "---0번째 Layer의 node 중에서 1번째 Layer의 3번째 노드와 엣지로 연결되어 있는 노드의 idx---" << endl;
    for (size_t i = 0; i < parent.size(); ++i) {
        cout << get<1>(parent[i]) << endl;
        // 해당 엣지 삭제
        directedGraph.removeEdge(parent[i], 3);
    }
    cout << "---위의 엣지를 제거한 후 graph 상태---" << endl;
    // 결과 확인용 
    directedGraph.printGraph();

    return 0;
}