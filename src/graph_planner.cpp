#include "graph_planner.hpp"
#include "NodeGraph.hpp"

void plotSpline(const Spline& spline, const string& color);

//////////////////////////////////////////////////////////////////////////////
// DVector dependent functions
//////////////////////////////////////////////////////////////////////////////

pair<DVector, DVector> computeBoundRight(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &width_r) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || width_r.empty()) {
        throw runtime_error("computeBoundRight() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_bound_r(len), y_bound_r(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_bound_r[i] = pos_x[i] + norm_x[i] * width_r[i];
        y_bound_r[i] = pos_y[i] + norm_y[i] * width_r[i];
    }
    
    return {x_bound_r, y_bound_r};

}

pair<DVector, DVector> computeBoundLeft(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &width_l) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || width_l.empty()) {
        throw runtime_error("computeBoundLeft() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_bound_l(len), y_bound_l(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_bound_l[i] = pos_x[i] - norm_x[i] * width_l[i];
        y_bound_l[i] = pos_y[i] - norm_y[i] * width_l[i];
    }
    
    return {x_bound_l, y_bound_l};

}

pair<DVector, DVector> computeRaceline(DVector &pos_x, DVector &pos_y,
                                         DVector &norm_x, DVector &norm_y,
                                         DVector &norm_l) {
    if (pos_x.empty() || pos_y.empty() || norm_x.empty() || norm_y.empty() || norm_l.empty()) {
        throw runtime_error("computeBoundRaceline() - Empty DVector !!");
    }

    int len = pos_x.size();
    DVector x_raceline(len), y_raceline(len);
    
    for (size_t i = 0; i < len; ++i) {
        x_raceline[i] = pos_x[i] + norm_x[i] * norm_l[i];
        y_raceline[i] = pos_y[i] + norm_y[i] * norm_l[i];
    }
    
    return {x_raceline, y_raceline};

}

DVector computeDeltaS(DVector &rl_s) {
    if (rl_s.empty()) {
        throw runtime_error("computeDeltaS() - Empty DVector !!");
    }

    int len = rl_s.size();
    DVector rl_ds(len);

    // 마지막 원소는 0
    for (size_t i = 0; i < len - 1; ++i) {
        rl_ds[i] = rl_s[i+1] - rl_s[i];
    }
    
    return rl_ds;

}

DVector computeHeading(DVector &x_raceline, DVector &y_raceline) {

    DVector psi;
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
    return psi;
}

DMap loadGlobalTrajectoryMap(string fname) {
  DMap gtMap = readDMapFromCSV(fname);

  auto [rb_x, rb_y] = computeBoundRight(gtMap[POS_X], gtMap[POS_Y],
                                        gtMap[NORM_X], gtMap[NORM_Y],
                                        gtMap[WIDTH_R]);
  gtMap[RB_X] = rb_x;
  gtMap[RB_Y] = rb_y;

  auto [lb_x, lb_y] = computeBoundLeft(gtMap[POS_X], gtMap[POS_Y],
                                       gtMap[NORM_X], gtMap[NORM_Y],
                                       gtMap[WIDTH_L]);
  gtMap[LB_X] = lb_x;
  gtMap[LB_Y] = lb_y;

  auto [rl_x, rl_y] = computeRaceline(gtMap[POS_X], gtMap[POS_Y],
                                      gtMap[NORM_X], gtMap[NORM_Y],
                                      gtMap[NORM_L]);
  gtMap[RL_X] = rl_x;
  gtMap[RL_Y] = rl_y;

  DVector rl_ds = computeDeltaS(gtMap[RL_S]);
  gtMap[RL_dS] = rl_ds;

  return gtMap;
}

IVector sampleLayersFromRaceline(const DVector& kappaVector,
                              const DVector& distVector,
                              const YAML::Node& params) {

    IVector layerIndexesSampled;
    const size_t n = kappaVector.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    // params
    float curve_thr = params["lattice"]["curve_thr"].as<float>();
    float d_straight = params["lattice"]["d_straight"].as<float>();
    float d_curve = params["lattice"]["d_curve"].as<float>();

    for (size_t i = 0; i < n; ++i) {
        // 곡선이면 최소 거리 갱신
        if ((cur_dist + distVector[i]) > next_dist_min && fabs(kappaVector[i]) > curve_thr) {
            next_dist = cur_dist;
        }
        // cout << fabs(kappaVector[i]) << endl;
        // 다음 샘플링 지점 도달
        if ((cur_dist + distVector[i]) > next_dist) {
            layerIndexesSampled.push_back(static_cast<int>(i));

            if (fabs(kappaVector[i]) < curve_thr) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += distVector[i];
    }

    // for (size_t i=0; i < layerIndexesSampled.size(); ++i) 
    //     cout << layerIndexesSampled[i] << endl;
    // cout << "size: " << layerIndexesSampled.size() << endl;

    return layerIndexesSampled;
}

DMap createSampledTrajectoryMap(DMap gtMap, YAML::Node params) {
  DMap stMap;
  IVector layerIndexesSampled = sampleLayersFromRaceline(gtMap[RL_KAPPA], gtMap[RL_dS], params);

  for (const auto& [key, vec] : gtMap) {
    for (int idx : layerIndexesSampled) {
      if (idx >= 0 && idx < vec.size()) {
        stMap[key].push_back(vec[idx]);
      } 
    }
  }

  stMap[RL_dS] = computeDeltaS(stMap[RL_S]);
  stMap[RL_PSI] = computeHeading(stMap[RL_X], stMap[RL_Y]);
  stMap[LB_PSI] = computeHeading(stMap[LB_X], stMap[LB_Y]);
  stMap[RB_PSI] = computeHeading(stMap[RB_X], stMap[RB_Y]);  

  return stMap;
}

auto createNodeMap(DMap &stMap, const YAML::Node &params) -> pair<NodeMap, IVector> {
    NodeMap nodeMap;
    IVector nodeIndexesOnRaceline;
    
    const int N = stMap[NORM_L].size();
    Vector2d node_pos;
    nodeMap.resize(N);    // N개 레이어 기준, nodeMap 벡터를 N 크기로 초기화 (각 레이어에 노드 저장)
    
    // params
    float veh_width = params["vehicle"]["veh_width"].as<float>();
    float lat_resolution = params["lattice"]["lat_resolution"].as<float>();
    
    // layer 별로 loop 돈다. for 루프 안이 한 레이어 내에서 하는 작업 내용물.
    for (size_t i = 0; i < N; ++i){ 
        Node node;
        // raceline이 layer 내에서 몇 번째 인덱스인지 확인. 이를 기준으로 node의 첫 번째 기준을 삼을 예정(s).
        int raceline_index = floor((stMap[WIDTH_L][i] + stMap[NORM_L][i] - veh_width / 2) / lat_resolution);
        nodeIndexesOnRaceline.push_back(raceline_index);

        Vector2d ref_xy(stMap[POS_X][i], stMap[POS_Y][i]);    // 기준선에서의 위치
        Vector2d norm_vec(stMap[NORM_X][i], stMap[NORM_Y][i]);  // 기준선에서 수직한 노멀 벡터 따라 노드 배치
        
        double start_alpha = stMap[NORM_L][i] - raceline_index * lat_resolution;    // 제일 왼쪽 노드가 노멀 벡터를 따라 얼마나 떨어져 있는지
        int node_idx = 0;
        int num_nodes = (stMap[WIDTH_R][i] + stMap[WIDTH_L][i] - veh_width) / lat_resolution;  // num_nodes : 좌우 총 가능한 노드 수
        
        nodeMap[i].resize(num_nodes); 
        
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

            nodeMap[i][node_idx] = node;
            ++node_idx;

        }

    }
    return {nodeMap, nodeIndexesOnRaceline};
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
    int closestIdx = -1;
    double min_dist2 = numeric_limits<double>::max();
    for (int i = 0; i < centerline.rows() - 1; ++i) {
        // segment 중심 계산
        Vector2d mid = (centerline.row(i) + centerline.row(i + 1)) / 2.0;
        double dist2 = (mid - pos).squaredNorm();
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closestIdx = i;
        }
    }

    if (closestIdx < 0 || closestIdx >= bound_l.rows() - 1)
        return false; // 예외 처리

    // bound_l, bound_r, centerline 보간 (선형 보간 10개 지점)
    int interp_points = 10;
    MatrixXd bl_interp(interp_points, 2);
    MatrixXd br_interp(interp_points, 2);
    MatrixXd center_interp(interp_points, 2);

    for (int i = 0; i < interp_points; ++i) {
        double t = static_cast<double>(i) / (interp_points - 1);
        bl_interp.row(i) = (1 - t) * bound_l.row(closestIdx) + t * bound_l.row(closestIdx + 1);
        br_interp.row(i) = (1 - t) * bound_r.row(closestIdx) + t * bound_r.row(closestIdx + 1);
        center_interp.row(i) = (1 - t) * centerline.row(closestIdx) + t * centerline.row(closestIdx + 1);
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

IPair getClosestNodes(const NodeMap& nodeMap, const Vector2d& pos, int limit=1) {
    IPair closestIdx;
    int num_nodes = 0;
    for (const auto& layer : nodeMap) {
        num_nodes += layer.size();
    }

    MatrixXd node_xy(num_nodes, 2);
    int idx = 0;

    for (size_t i = 0; i < nodeMap.size(); ++i) {
        for (size_t j = 0; j < nodeMap[i].size(); ++j) {
            const Node& node = nodeMap[i][j];
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
    for (size_t i = 0; i < nodeMap.size(); ++i) {
        for (size_t j = 0; j < nodeMap[i].size(); ++j) {
            dist_info.emplace_back(dist2(re_idx++), i, j);
        }
    }

    // 최소 거리 limit개만 앞으로 정렬
    nth_element(dist_info.begin(), dist_info.begin() + limit, dist_info.end());

    // 결과 저장
    for (int k = 0; k < limit; ++k) {
        auto [dist, i, j] = dist_info[k];
        closestIdx = make_pair(i, j);
        cout << "Closest node: layer=" << i << ", idx=" << j << endl;
    }
    return closestIdx;
}

auto computeSplines(const MatrixXd &path,
                    double psi_s,
                    double psi_e,
                    bool use_dist_scaling) -> unique_ptr<Spline> {
    // 구간 길이 계산
    VectorXd el_lengths;
    if (use_dist_scaling)
    {
        el_lengths.resize(path.rows() - 1);
        for (int i = 0; i < path.rows() - 1; ++i)
        {
            el_lengths(i) = (path.row(i + 1) - path.row(i)).norm();
        }
    }
    // 맨 마지막 거리 추가
    if (use_dist_scaling)
    {
        VectorXd el_tmp(el_lengths.size() + 1);
        el_tmp << el_lengths, el_lengths(0);
        el_lengths = el_tmp;
    }

    // 도함수 스케일링
    // 인접 구간 간 거리 비율로 스케일링 계수를 만들어 도함수 연속 조건 맞춤
    int no_splines = path.rows() - 1;
    VectorXd scaling = VectorXd::Ones(no_splines - 1);
    if (use_dist_scaling)
    {
        for (int i = 0; i < no_splines - 1; ++i)
        {
            scaling(i) = el_lengths(i) / el_lengths(i + 1);
        }
    }

    MatrixXd M = MatrixXd::Zero(no_splines * 4, no_splines * 4);
    VectorXd b_x = VectorXd::Zero(no_splines * 4);
    VectorXd b_y = VectorXd::Zero(no_splines * 4);

    // spline 위치/도함수/2차 도함수 연속 조건 표현
    Matrix<double, 4, 8> template_M;
    template_M << 1, 0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 0, 0, 0, 0,
        0, 1, 2, 3, 0, -1, 0, 0,
        0, 0, 2, 6, 0, 0, -2, 0;
    // spline 구간별로 행렬 세팅
    // 마지막 spline은 위치조건만
    for (int i = 0; i < no_splines; ++i)
    {
        int j = i * 4;
        if (i < no_splines - 1)
        {
            M.block(j, j, 4, 8) = template_M;
            M(j + 2, j + 5) *= scaling(i);         // 이웃한 구간 도함수 값 일치하도록
            M(j + 3, j + 6) *= pow(scaling(i), 2); // 2차 도함수 가속도 연속 조건
        }
        else
        {
            M.block(j, j, 2, 4) << 1, 0, 0, 0,
                1, 1, 1, 1;
        }
        // b_x.segment(j, 2): b_x[j]와 b_x[j+1]
        b_x.segment(j, 2) << path(i, 0), path(i + 1, 0); // x좌표에 대한 위치 조건 벡터
        b_y.segment(j, 2) << path(i, 1), path(i + 1, 1); // y좌표에 대한 위치 조건 벡터
    }
    // 시작/끝점에서의 psi 반영
    psi_s += M_PI_2;
    psi_e += M_PI_2;

    M(no_splines * 4 - 2, 1) = 1.0;
    double el_length_s = el_lengths.size() > 0 ? el_lengths(0) : 1.0;
    b_x(no_splines * 4 - 2) = cos(psi_s) * el_length_s;
    b_y(no_splines * 4 - 2) = sin(psi_s) * el_length_s;
    // el_lengths.tail(1): 벡터의 마지막 원소 반환
    // el_lengths.tail(1)(0): 원소 값 가져옴
    M.block(no_splines * 4 - 1, no_splines * 4 - 4, 1, 4) << 0, 1, 2, 3;
    double el_length_e = el_lengths.size() > 0 ? el_lengths.tail(1)(0) : 1.0;
    // 끝점에서의 곡선의 방향(psi_e)이 실제 경로의 마지막 구간 길이에 맞게 변화량이 되도록 스케일을 맞춰주는 것
    // 방향만 넣으면, 구간 길이가 1로 가정된 것처럼 도함수 조건이 설정되어 실제 경로의 스케일과 맞지 않게 된다.
    b_x(no_splines * 4 - 1) = cos(psi_e) * el_length_e;
    b_y(no_splines * 4 - 1) = sin(psi_e) * el_length_e;

    VectorXd x_les = M.fullPivLu().solve(b_x);
    VectorXd y_les = M.fullPivLu().solve(b_y);
    MatrixXd coeffs_x = x_les.transpose();
    MatrixXd coeffs_y = y_les.transpose();
    // cout << coeffs_x.cols() << endl;

    VectorXd kappaVector;
    vector<Vector2d> points_xy;
    double cost = 0.0;
    bool raceline = false;
    // 결과 반환
    return make_unique<Spline>(Spline{
        coeffs_x, // (4, 1)
        coeffs_y, // (4, 1)
        kappaVector,
        el_lengths,
        points_xy,
        cost,
        raceline,
    });
}

auto sampleSingleSpline(MatrixXd &coeffs_x, MatrixXd &coeffs_y, YAML::Node &params) -> pair<vector<Vector2d>, VectorXd> {

    if (coeffs_x.rows() != coeffs_y.rows())
    {
        throw invalid_argument("Coefficient matrices must have the same length!");
    }

    if (coeffs_x.cols() == 2 && coeffs_y.cols() == 2)
    {
        throw invalid_argument("Coefficient matrices do not have two dimensions!");
    }
    int no_interp_points = params["sampling"]["no_interp_points"].as<int>();

    VectorXd t_steps(no_interp_points);
    double step = 1.0 / (no_interp_points - 1);
    for (size_t i = 0; i < no_interp_points; ++i)
    {
        t_steps[i] = i * step;
    }

    vector<Vector2d> points_xy;
    VectorXd kappaVector(no_interp_points + 1);
    // kappaVector.reserve(no_interp_points + 1);

    for (int i = 0; i < no_interp_points; ++i)
    {
        double t = t_steps(i);
        double t2 = t * t;
        double t3 = t2 * t;

        // 좌표 계산
        double x = coeffs_x(0, 0) + coeffs_x(0, 1) * t + coeffs_x(0, 2) * t2 + coeffs_x(0, 3) * t3;
        double y = coeffs_y(0, 0) + coeffs_y(0, 1) * t + coeffs_y(0, 2) * t2 + coeffs_y(0, 3) * t3;

        // 1차 미분
        double x_d = coeffs_x(0, 1) + 2 * coeffs_x(0, 2) * t + 3 * coeffs_x(0, 3) * t2;
        double y_d = coeffs_y(0, 1) + 2 * coeffs_y(0, 2) * t + 3 * coeffs_y(0, 3) * t2;

        // 2차 미분
        double x_dd = 2 * coeffs_x(0, 2) + 6 * coeffs_x(0, 3) * t;
        double y_dd = 2 * coeffs_y(0, 2) + 6 * coeffs_y(0, 3) * t;

        double denom = pow(x_d * x_d + y_d * y_d, 1.5);

        kappaVector(i) = (x_d * y_dd - y_d * x_dd) / denom;
        points_xy.emplace_back(x, y);
    }

    return {points_xy, kappaVector};
}

void setInitialPose(DMap &stMap,
                    const NodeMap &nodeMap,
                    const IVector &nodeIndexesOnRaceline,
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
    
    IPair closestIdx = getClosestNodes(nodeMap, initial_pos);
    int start_layer = closestIdx.first;
    int start_node = closestIdx.second;
    double start_heading = nodeMap[closestIdx.first][closestIdx.second].psi;

    int end_layer = (closestIdx.first + initial_layer) % (nodeMap.size() - 1);

    for (int layer_idx = start_layer; layer_idx < end_layer;++layer_idx) {
        if (layer_idx != start_layer) {
            start_node = nodeIndexesOnRaceline[layer_idx];
            start_heading = nodeMap[layer_idx][start_node].psi;
        }

        int goal_layer = layer_idx + 1;
        int goal_node = nodeIndexesOnRaceline[goal_layer];
        double goal_heading = nodeMap[goal_layer][goal_node].psi;
        double heading_diff = abs(start_heading - goal_heading);
        Vector2d start_pos(nodeMap[layer_idx][start_node].x, nodeMap[layer_idx][start_node].y);
        Vector2d end_pos(nodeMap[goal_layer][goal_node].x, nodeMap[goal_layer][goal_node].y);

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

        auto result = computeSplines(path, start_heading, goal_heading); // coeffs_x, coeffs_y, kappa, el_lengths, cost
        if (!result) {
            cerr << "computeSplines failed for edge: " << endl;
        continue;
    }

        auto [points_xy, kappaVector] = sampleSingleSpline(result->coeffs_x, result->coeffs_y, params);

        if (kappaVector.size() == 0 || points_xy.size() == 0) {
            throw invalid_argument("Points or kappaVector's Size is zero!! - interpSplines()");
        }
        result->kappaVector = kappaVector;
        result->points_xy = points_xy;

        plotSpline(*result, "blue");    
    }
    cout << "Goal node: layer=" << end_layer << ", idx= " << nodeIndexesOnRaceline[end_layer] << endl;
    
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