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
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
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
            IVector& raceline_index_array,
            const double veh_width,
            float lat_resolution) {
    
    const size_t N = sampling_map[__alpha].size();
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

// -------------------------------------------
// gen_spline --------------------------------
// -------------------------------------------


// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;          // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;          // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                 // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;// 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

// 경로 행렬의 각 점 사이의 유클리드 거리 계산 함수
VectorXd computeEuclideanDistances(const MatrixXd& path) {
    int N = path.rows() - 1;    // 구간 개수는 점 개수 - 1
    VectorXd dists(N);          // 거리 저장 벡터
    for (int i = 0; i < N; ++i) {
        dists(i) = (path.row(i + 1) - path.row(i)).norm();  // 두 점 사이 거리 계산 (유클리드 노름)
    }
    return dists;
}

// path : x,y 좌표 배열
// el_lengths 전달 안 하면 nullptr (점 사이 거리는 필수 x)
// psi_s, psi_e도 전달 안 하면 NAN (필수 x)
SplineResult calcSplines(const MatrixXd& path, const VectorXd* el_lengths = nullptr,
                         double psi_s = NAN, double psi_e = NAN,
                         bool use_dist_scaling = true) {

        // 1.폐곡선 여부 판단
        MatrixXd updated_path = path;
        
        // 시작 점과 끝 점 사이 거리가 매우 작으면 같다고 간주한다.
        bool closed = (path.row(0) - path.row(path.rows() - 1)).norm() < 1e-6;

        if (closed) {
            // 폐곡선 유지 : 시작점 다시 한 번 넣기 (마지막 스플라인 조각을 만들기 위해서)
            updated_path.conservativeResize(path.rows() + 1, Eigen::NoChange);
            updated_path.row(updated_path.rows() - 1) = path.row(0);
        }

        int no_splines = updated_path.rows() - 1;   // 경로 구간 수 = 점 개수 - 1

        VectorXd ds;
        // 유클리드 거리 담고 있는 el_lengths (nullptr -> 유클리드 거리 계산, 아니면 el_lengths가 가리키는 VerctorXd 객체 거리 사용)
        if (el_lengths == nullptr) {
            ds = computeEuclideanDistances(updated_path);  // 기본 거리 계산 함수 호출
        } else {
            ds = *el_lengths;  // 외부에서 받은 길이 벡터 사용
        }

        // 유클리드 거리 출력
        cout << "ds(거리벡터):\n" << ds.transpose() << "\n";

        // 2. 행렬 M 구성
        // 각 스플라인은 x, y 각각에 대해 4개의 계수(a0 ~ a3) 가짐
        // 스플라인 0 -> 행 0 ~ 3, 스플라인 1 -> 행 4 ~ 7

        // 위치 조건(각 스플라인은 시작점과 끝점에서 위치가 정확히 맞아야 함) : 식 2*(N-1)개
        // 1차 도함수 연속(기울기 연속) : 식 (N-2)개
        // 2차 도함수 연속(곡률 연속) : 식 (N-2)개
        // 양 끝 점(경계 조건) : 식 2개

        // 스케일링 없이 진행. 1, 2차 도합수 조건에 직접 거리 d 사용하면 됨
        
        MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
        VectorXd b_x = VectorXd::Zero(no_splines * 4);
        VectorXd b_y = VectorXd::Zero(no_splines * 4);

        // 위치 조건
        for (int i = 0; i < no_splines; i++) {
            int row_s = i * 4;  // 스플라인 시작 행
            double d = ds(i);

            // 시작점 위치 조건
            M(row_s, row_s) = 1.0;  // a_0 = p_i
            b_x(row_s) = updated_path(i,0);
            b_y(row_s) = updated_path(i,1);

            // 끝점 위치 조건
            RowVector4d T(1.0, d, d*d, d*d*d);
            M.block(row_s + 1, row_s, 1, 4) = T; // 시작 행, 시작 열, 선택할 행, 선택할 열
            b_x(row_s + 1) = updated_path(i+1, 0);
            b_y(row_s + 1) = updated_path(i+1, 1);
        }

        // 기울기 조건
        // i = 0 ~ no_splines - 2까지만 적용 (마지막 스플라인 뒤에는 비교할 스플라인이 없기 때문)
        for (int i = 0; i < no_splines - 1; ++i) {
            int row = 2 * no_splines + i;
            double d = ds(i);

            M(row, 4 * i + 1) = 1.0;              // a_1^i
            M(row, 4 * i + 2) = 2.0 * d;          // 2a_2^i * d
            M(row, 4 * i + 3) = 3.0 * d * d;      // 3a_3^i * d^2
            M(row, 4 * (i + 1) + 1) = -1.0;       // -a_1^{i+1}

            b_x(row) = 0.0;
            b_y(row) = 0.0;
        }

        // 곡률 조건
        for (int i = 0; i < no_splines - 1; ++i) {
            int row = 3 * no_splines - 1 + i;
            double d = ds(i);

            M(row, 4 * i + 2) = 2.0;                // 2a_2^i
            M(row, 4 * i + 3) = 6.0 * d;            // 6a_3^i * d
            M(row, 4 * (i + 1) + 2) = -2.0;         // -2a_2^{i+1}

            b_x(row) = 0.0;
            b_y(row) = 0.0;
        }

        // 닫히지 않은 경로에서 시작점과 끝점 방향 조건
        if (!closed) {
            // 시작점 도함수 조건 (t=0)
            int start_row = 4 * no_splines - 2; // 마지막에서 두 번째 줄이 시작점 도함수 조건
            M.row(start_row).setZero();
            M(start_row, 1) = 1.0;  // a_1항에 1 설정 (S'(0) = a_1)

            double d_s = ds(0); // 첫 번째 구간 길이
            b_x(start_row) = cos(psi_s) * d_s;  // x방향 단위벡터 * d_s
            b_y(start_row) = sin(psi_s) * d_s;  // y방향 단위벡터 * d_s

            // 끝점 도함수 조건 (t=d)
            int end_row = 4 * no_splines - 1;
            M.row(end_row).setZero();
            
            double d_e = ds(no_splines - 1);    // 마지막 구간 길이
            int offset = 4 * (no_splines - 1);  // 마지막 스플라인 계수들의 시작 인덱스
            M(end_row, offset + 1) = 1.0;
            M(end_row, offset + 2) = 2.0 * d_e;
            M(end_row, offset + 3) = 3.0 * d_e * d_e;

            b_x(end_row) = cos(psi_e) * d_e;
            b_y(end_row) = sin(psi_e) * d_e;
        }

        // [DEBUG] M, b_x, b_y 출력
        cout << "[DEBUG] M 행렬:\n" << M << "\n";
        cout << "[DEBUG] b_x 벡터:\n" << b_x.transpose() << "\n";
        cout << "[DEBUG] b_y 벡터:\n" << b_y.transpose() << "\n";

        // 연립방정식 풀기
        VectorXd x_les = M.colPivHouseholderQr().solve(b_x);
        VectorXd y_les = M.colPivHouseholderQr().solve(b_y);

        MatrixXd coeffs_x = Map<MatrixXd>(x_les.data(), 4, no_splines).transpose();
        MatrixXd coeffs_y = Map<MatrixXd>(y_les.data(), 4, no_splines).transpose();


        // 법선 벡터 계산 (노멀 벡터)
        MatrixXd normvec(no_splines, 2);
        for (int i = 0; i < no_splines; ++i) {
            double dx = coeffs_x(i, 1);  // a₁: x방향 도함수 시작점
            double dy = coeffs_y(i, 1);  // a₁: y방향 도함수 시작점
            normvec(i, 0) = -dy;
            normvec(i, 1) = dx;
        }

       // 정규화된 법선 벡터 계산
        MatrixXd normvec_normalized(no_splines, 2);
        for (int i = 0; i < no_splines; ++i) {
            double norm = normvec.row(i).norm();
            if (norm > 1e-6)
                normvec_normalized.row(i) = normvec.row(i) / norm;
            else
                normvec_normalized.row(i) = Vector2d(0.0, 0.0); // norm이 너무 작을 때 0으로 대체
        }

        // 결과 구조체에 저장
        SplineResult result;
        result.coeffs_x = coeffs_x;
        result.coeffs_y = coeffs_y;
        result.M = M;
        result.normvec_normalized = normvec_normalized;


        return result;

}

Vector2d evaluateSpline(const RowVector4d& coeffs_x, const RowVector4d& coeffs_y, double t) {
    double x = coeffs_x(0) + coeffs_x(1) * t + coeffs_x(2) * t * t + coeffs_x(3) * t * t * t;
    double y = coeffs_y(0) + coeffs_y(1) * t + coeffs_y(2) * t * t + coeffs_y(3) * t * t * t;
    return Vector2d(x, y);
}

void visualizeSplines(const SplineResult& result, const MatrixXd& path) {
    using namespace matplotlibcpp;

    vector<double> xs, ys;
    int samples_per_segment = 20;

    for (int i = 0; i < result.coeffs_x.rows(); ++i) {
        double d = (path.row(i + 1) - path.row(i)).norm();

        for (int j = 0; j <= samples_per_segment; ++j) {
            double t = d * j / static_cast<double>(samples_per_segment);
            Vector2d pt = evaluateSpline(result.coeffs_x.row(i), result.coeffs_y.row(i), t);
            xs.push_back(pt.x());
            ys.push_back(pt.y());
        }
    }

    // 원래 입력 경로 점
    vector<double> x_orig, y_orig;
    for (int i = 0; i < path.rows(); ++i) {
        x_orig.push_back(path(i, 0));
        y_orig.push_back(path(i, 1));
    }

    figure();
    plot(xs, ys, "r-");           // 스타일만 사용
    plot(x_orig, y_orig, "bo--"); // 스타일만 사용
    legend();                     // 레전드 호출 시 자동 생성
    title("Spline Visualization");
    axis("equal");
    grid(true);
    show();
}



void genEdges() {
    
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

    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);


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

    SplineResult result = calcSplines(path_xy, nullptr, psi_s, psi_e);

    // 시각화
    visualizeSplines(result, path_xy);


    //genEdges();
    // visual process 
    //visual();

    return 0;
}
