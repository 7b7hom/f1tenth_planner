#include "graph_planner.hpp"

unique_ptr<SplineResult> calcSplines(const MatrixXd &path,
                                     double psi_s = NAN,
                                     double psi_e = NAN,
                                     bool use_dist_scaling = true) {
    // 구간 길이 계산
    VectorXd el_lengths;
    if (use_dist_scaling) {
        el_lengths.resize(path.rows() - 1);
        for (int i = 0; i < path.rows() - 1; ++i) {
            el_lengths(i) = (path.row(i+1) - path.row(i)).norm();
        }
    } 
    // 맨 마지막 거리 추가
    if (use_dist_scaling) {
        VectorXd el_tmp(el_lengths.size() + 1);
        el_tmp << el_lengths, el_lengths(0);
        el_lengths = el_tmp;
    }
    // 도함수 스케일링
    // 인접 구간 간 거리 비율로 스케일링 계수를 만들어 도함수 연속 조건 맞춤
    int no_splines = path.rows() - 1;
    VectorXd scaling = VectorXd::Ones(no_splines - 1);
    if (use_dist_scaling) {
        for (int i = 0; i < no_splines - 1; ++i) {
            scaling(i) = el_lengths(i) / el_lengths(i+1);
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
    for (int i = 0; i < no_splines; ++i) {
        int j = i * 4;
        if (i < no_splines - 1) {
            M.block(j, j, 4, 8) = template_M;
            M(j+2, j+5) *= scaling(i); // 이웃한 구간 도함수 값 일치하도록 
            M(j+3, j+6) *= pow(scaling(i), 2); // 2차 도함수 가속도 연속 조건
        } else {
            M.block(j, j, 2, 4) << 1, 0, 0, 0,
                                   1, 1, 1, 1;
        }
        // b_x.segment(j, 2): b_x[j]와 b_x[j+1]
        b_x.segment(j, 2) << path(i, 0), path(i+1, 0); // x좌표에 대한 위치 조건 벡터
        b_y.segment(j, 2) << path(i, 1), path(i+1, 1); // y좌표에 대한 위치 조건 벡터
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

    // 결과 반환
    return make_unique<SplineResult>(SplineResult{
        coeffs_x,  // (4, 1)
        coeffs_y,  // (4, 1)
        el_lengths,
    });
}
// 단일 스플라인에 대한 spline에 대한 샘플링 
MatrixXd* interpSplines(MatrixXd &coeffs_x,
                   MatrixXd &coeffs_y,
                   float stepsize_approx,
                   double spline_len = NAN, 
                   int no_interp_points = 15) {
    if (coeffs_x.rows() != coeffs_y.rows()) {
        throw invalid_argument("Coefficient matrices must have the same length!");
    }

    if (coeffs_x.cols() == 2 && coeffs_y.cols() == 2) {
        throw invalid_argument("Coefficient matrices do not have two dimensions!");
    }

    if (isnan(stepsize_approx)) {
        throw invalid_argument("Provide one of 'stepsize_approx' and 'stepnum_fixed' and set the other to 'None'!");
    }
    // spline 계산 
    // 샘플링할 점 개수
    if (isnan(spline_len)) {
        int no_splines = coeffs_x.rows();
        // cout << no_splines << endl;
        VectorXd t_steps(no_interp_points);
        double step = 1.0 / (no_interp_points - 1);

        for (size_t i = 0; i < no_interp_points; ++i) {
            t_steps[i] = i*step;
        }

        MatrixXd* spl_coords = new MatrixXd(no_interp_points, 2);

        for (int i = 0; i < no_splines; ++i) {
            spl_coords->col(0) =
                coeffs_x(i, 0) * VectorXd::Ones(no_interp_points)
                + coeffs_x(i, 1) * t_steps
                + coeffs_x(i, 2) * t_steps.array().pow(2).matrix()
                + coeffs_x(i, 3) * t_steps.array().pow(3).matrix();

            spl_coords->col(1) =
                coeffs_y(i, 0) * VectorXd::Ones(no_interp_points)
                + coeffs_y(i, 1) * t_steps
                + coeffs_y(i, 2) * t_steps.array().pow(2).matrix()
                + coeffs_y(i, 3) * t_steps.array().pow(3).matrix();           
        }
        spline_len = 0.0;
        for (int j = 1; j < no_interp_points; ++j) {
            double dx = (*spl_coords)(j, 0) - (*spl_coords)(j-1, 0);
            double dy = (*spl_coords)(j, 1) - (*spl_coords)(j-1, 1);

            spline_len +=sqrt(dx*dx + dy*dy);
        }
        // cout << spline_len << endl; // levine: 0.5~2.5
        return spl_coords;
    }

    return nullptr;
    
    
    // 보류
    // MatrixXd path_interp = MatrixXd::Zero(no_interp_points, 2); 
    // VectorXd spline_inds = VectorXd::Zero(no_interp_points);
    // VectorXd t_values = VectorXd::Zero(no_interp_points);

    // for (size_t i = 0; i < no_interp_points - 1; ++i) {
    // }
    
}

#if 0
bool isInsideTrack() {

}
#endif

void genEdges(NodeMap &nodesPerLayer, 
              Graph &edgeList,
              SplineMap &splineMap,
              const IVector &raceline_index_array,
              const float lat_offset,
              const float lat_resolution,
              const float curve_thr,
              const int max_lat_steps,
              const float stepsize_approx) {
    
    if (lat_offset <= 0.0) {
        throw invalid_argument("Too small lateral offset!");
    }

    // cout << nodesPerLayer.size() << endl; 출력: 51
    // 레이어 별 loop
    for (int layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) {
        
        int start_layer = layerIdx;
        int end_layer = layerIdx + 1;

        // cout << "start_layer:" << start_layer << endl;
        // cout << "nodesPerLayer.size()" << nodesPerLayer.size() << endl;

        // 마지막 layer의 경우 0번째 layer와 연결시킬 수 있도록 end_layer 조정 
        if (end_layer >= nodesPerLayer.size()) {
            end_layer -= nodesPerLayer.size();
        }

        // start layer 내 노드별 loop
        for (size_t startIdx = 0; startIdx < nodesPerLayer[start_layer].size(); ++startIdx) {
            // 기준 노드
            Node &startNode = nodesPerLayer[start_layer][startIdx];
            
            int refDestIdx = raceline_index_array[end_layer] - (raceline_index_array[start_layer] - startIdx);
            refDestIdx = max(0, min(refDestIdx, static_cast<int>(nodesPerLayer[end_layer].size() -1)));
            // refDestIdx = clamp(refDestIdx, 0, static_cast<int>(nodesPerLayer[end_layer].size() - 1));
            // int refDestIdx = startIdx;
            // int refDestIdx = raceline_index_array[end_layer];
            Node &srcEndNode = nodesPerLayer[end_layer][refDestIdx];

            Vector2d d_start(startNode.x, startNode.y);
            Vector2d d_end(srcEndNode.x, srcEndNode.y);

            // spline 연결할 노드 선정 기준 : lat_steps
            double dist = (d_end - d_start).norm();
            // genNode에서 kappa 계산한거 토대로(+기능 추가 완료)
            double ratio = min(startNode.kappa / curve_thr, 2.0); // 최대 2배까지만 증폭
            double factor = 1.0 + 0.5 * ratio;
            // 커브에서 더 많이 연결(추월 경로를 위하여)
            int lat_steps = round(factor * dist * lat_offset / lat_resolution);
            // cout << start_layer << "의 " << startIdx << "가 다음 refendNode와의 거리: " << dist << endl;
            lat_steps = min(lat_steps, max_lat_steps); // srcEndNode 기준 2*lat_steps + 1개의 노드와 연결한다.
            // cout << startIdx << "번째 노드의 lat_steps" << lat_steps << endl;
            // startNode와 lat_steps 기준 해당되는 노드들 spline 연결 
            for (int destIdx = max(0, refDestIdx - lat_steps); 
                destIdx <= min(static_cast<int>(nodesPerLayer[end_layer].size() - 1), refDestIdx + lat_steps); ++destIdx) {
                    Node &endNode = nodesPerLayer[end_layer][destIdx];
                    
                    MatrixXd path(2, 2);
                    path(0,0) = startNode.x;
                    path(0,1) = startNode.y;
                    path(1,0) = endNode.x;
                    path(1,1) = endNode.y;

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
    int edge_cnt = 0;
    // layer 개수만큼 loop
    for (size_t i = 0; i < nodesPerLayer.size();++i) {
      int start_layer = i;
      // start_layer에서의 노드 개수만큼 loop
      for (size_t s = 0; s < nodesPerLayer[start_layer].size(); ++s) {
        IPairVector childNode;
        IPair start = make_pair(i ,s);
        edgeList.getChildIdx(start, childNode);
        for (const auto& child : childNode) {
          edge_cnt++;
          EdgeKey srcKey = make_pair(start, child);
          // int end_layer = child.first;
          // int e = child.second;

          MatrixXd coeffs_x = splineMap[srcKey].coeffs_x;
          MatrixXd coeffs_y = splineMap[srcKey].coeffs_y;

          MatrixXd* spl_coords = interpSplines(coeffs_x, coeffs_y, stepsize_approx);
        

          delete spl_coords;
        //   calcHeadCurv(coeffs_x, coeffs_y);
        }
      }
    }
}

