#include "graph.hpp"
// #include "graph_planner.hpp"
class SplineHandler {
private:
    SplineMap splineMap;

public:
    SplineMap& getSplineMap() { return splineMap; }

    Spline& at(const IPair& start, const IPair& end) {
        return splineMap[start][end];
    }

    const Spline& at(const IPair& start, const IPair& end) const {
        return splineMap.at(start).at(end);
    }

    void writeSplineMapToCSV(const string& filename) {
        ofstream fout(filename);
        if (!fout.is_open()) throw runtime_error("Cannot open file");

        // 헤더
        fout << "start_layer,start_idx,end_layer,end_idx,coeffs_x,coeffs_y,kappa,points_x,points_y,raceline\n";

        for (const auto& [start, endMap] : splineMap) {
            for (const auto& [end, spline] : endMap) {
                fout << start.first << "," << start.second << ","
                    << end.first << "," << end.second << ",";

                // coeffs_x
                for (int i = 0; i < spline.coeffs_x.rows(); ++i)
                    for (int j = 0; j < spline.coeffs_x.cols(); ++j)
                        fout << spline.coeffs_x(i,j) << (i==spline.coeffs_x.rows()-1 && j==spline.coeffs_x.cols()-1 ? "," : " ");

                // coeffs_y
                for (int i = 0; i < spline.coeffs_y.rows(); ++i)
                    for (int j = 0; j < spline.coeffs_y.cols(); ++j)
                        fout << spline.coeffs_y(i,j) << (i==spline.coeffs_y.rows()-1 && j==spline.coeffs_y.cols()-1 ? "," : " ");

                // kappa
                for (int i = 0; i < spline.kappa.size(); ++i)
                    fout << spline.kappa(i) << (i==spline.kappa.size()-1 ? "," : " ");

                // points_xy
                for (int i = 0; i < spline.points_xy.size(); ++i)
                    fout << spline.points_xy[i].x() << " " << spline.points_xy[i].y() 
                        << (i==spline.points_xy.size()-1 ? "," : " ");

                fout << (spline.raceline ? 1 : 0) << "\n";
            }
        }
        fout.close();
    }

    void readSplineMapFromCSV(const string& filename) {
        ifstream fin(filename);
        if (!fin.is_open()) throw runtime_error("Cannot open file");

        string line;
        getline(fin, line); // 헤더 스킵

        while (getline(fin, line)) {
            stringstream ss(line);
            string item;

            IPair start, end;
            Spline spline;

            // start_layer, start_idx, end_layer, end_idx
            getline(ss, item, ','); start.first = stoi(item);
            getline(ss, item, ','); start.second = stoi(item);
            getline(ss, item, ','); end.first = stoi(item);
            getline(ss, item, ','); end.second = stoi(item);

            // coeffs_x
            getline(ss, item, ',');
            stringstream sx(item);
            DVector vx;
            double val;
            while (sx >> val) vx.push_back(val);

            spline.coeffs_x = MatrixXd(4, 1);
            for (int i=0;i<4;i++) for(int j=0;j<1;j++) spline.coeffs_x(i,j)=vx[i*1+j];

            // coeffs_y
            getline(ss, item, ',');
            stringstream sy(item);
            DVector vy;
            while (sy >> val) vy.push_back(val);
            spline.coeffs_y = MatrixXd(4 ,1);
            for (int i=0;i<4;i++) for(int j=0;j<1;j++) spline.coeffs_y(i,j)=vy[i*1+j];

            // kappa
            getline(ss, item, ',');
            stringstream sk(item);
            DVector vk;
            while (sk >> val) vk.push_back(val);
            spline.kappa = VectorXd::Map(vk.data(), vk.size());

            // points_xy
            getline(ss, item, ',');
            stringstream sp(item);
            vector<Vector2d> pts;
            double x,y;
            while (sp >> x >> y) pts.emplace_back(x,y);
            spline.points_xy = pts;

            // raceline
            getline(ss, item, ',');
            spline.raceline = (stoi(item) != 0);

            splineMap[start][end] = spline;
        }

        fin.close();
    }

    auto calcSplines(MatrixXd &path,
                     double psi_s,
                     double psi_e,
                     bool use_dist_scaling = true) -> unique_ptr<Spline> {
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

        VectorXd kappa;
        vector<Vector2d> points_xy;
        double cost = 0.0;
        bool raceline = false;
        // 결과 반환
        return make_unique<Spline>(Spline{
            coeffs_x,  // (4, 1)
            coeffs_y,  // (4, 1)
            kappa,
            el_lengths,
            points_xy,
            cost,
            raceline,
        });
    }

    auto genEdges(DMap &stMap,
                NodeMap &nodesPerLayer,
                const IVector &raceline_index_array,
                YAML::Node &params) -> Graph {

        Graph wayptGraph; // a graph of waypoints

        float lat_offset = params["lattice"]["lat_offset"].as<float>();
        float lat_resolution = params["lattice"]["lat_resolution"].as<float>();
        int max_lat_steps = params["lattice"]["max_lat_steps"].as<int>();

        if (lat_offset <= 0.0) {
            throw invalid_argument("Too small lateral offset!");
        }

        // raceline 
        for (int layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) { 
            int dstLayerIdx = layerIdx+1;
            if (dstLayerIdx >= nodesPerLayer.size()) {
                dstLayerIdx -= nodesPerLayer.size();
            }
            int startNodeIdx = raceline_index_array[layerIdx];
            int endNodeIdx = raceline_index_array[dstLayerIdx];

            Node& startNode = nodesPerLayer[layerIdx][startNodeIdx];
            Node& endNode = nodesPerLayer[dstLayerIdx][endNodeIdx];

            MatrixXd path(2, 2);
            path(0,0) = startNode.x;
            path(0,1) = startNode.y;
            path(1,0) = endNode.x;
            path(1,1) = endNode.y;

            auto result = this->calcSplines(path, startNode.psi, endNode.psi);

            IPair startPoint = make_pair(layerIdx, startNodeIdx);
            IPair endPoint = make_pair(dstLayerIdx, endNodeIdx);
            
            result->raceline = true;

            splineMap[startPoint][endPoint] = *result;

            wayptGraph.addEdge(startPoint, endPoint);
        }

        // cout << nodesPerLayer.size() << endl; 출력: 51
        // 레이어 별 loop
        // raceline spline 먼저 생성해서 splineMap에 등록, 뒤에서 재등록하지 않게 index 겹치면 pass
        for (int layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) {
            
            int srcLayerIdx = layerIdx;
            int dstLayerIdx = layerIdx + 1;

            // cout << "srcLayerIdx:" << srcLayerIdx << endl;
            // cout << "nodesPerLayer.size()" << nodesPerLayer.size() << endl;

            // 마지막 layer의 경우 0번째 layer와 연결시킬 수 있도록 dstLayerIdx 조정 
            if (dstLayerIdx >= nodesPerLayer.size()) {
                dstLayerIdx -= nodesPerLayer.size();
            }

            // start layer 내 노드별 loop
            for (size_t srcNodeIdx = 0; srcNodeIdx < nodesPerLayer[srcLayerIdx].size(); ++srcNodeIdx) {
                // 기준 노드
                Node& startNode = nodesPerLayer[srcLayerIdx][srcNodeIdx];
                
                int refEndNodeIdx = raceline_index_array[dstLayerIdx] - (raceline_index_array[srcLayerIdx] - srcNodeIdx);
                refEndNodeIdx = max(0, min(refEndNodeIdx, static_cast<int>(nodesPerLayer[dstLayerIdx].size() -1)));
                // refEndNodeIdx = clamp(refEndNodeIdx, 0, static_cast<int>(nodesPerLayer[dstLayerIdx].size() - 1));
                // int refEndNodeIdx = srcNodeIdx;
                // int refEndNodeIdx = raceline_index_array[dstLayerIdx];
                Node& endNode = nodesPerLayer[dstLayerIdx][refEndNodeIdx];

                Vector2d d_start(startNode.x, startNode.y);
                Vector2d d_end(endNode.x, endNode.y);

                // spline 연결할 노드 선정 기준 : lat_steps
                double dist = (d_end - d_start).norm();
                // genNode에서 kappa 계산한거 토대로(+기능 추가 완료)

                int lat_steps = static_cast<int>(round(dist * lat_offset / lat_resolution));
                // cout << srcLayerIdx << "의 " << srcNodeIdx << "가 다음 refendNode와의 거리: " << dist << endl;
                lat_steps = min(lat_steps, max_lat_steps); // endNode 기준 2*lat_steps + 1개의 노드와 연결한다.
                // cout << srcNodeIdx << "번째 노드의 lat_steps" << lat_steps << endl;
                // startNode와 lat_steps 기준 해당되는 노드들 spline 연결 
                for (int endNodeIdx = max(0, refEndNodeIdx - lat_steps); 
                    endNodeIdx <= min(static_cast<int>(nodesPerLayer[dstLayerIdx].size() - 1), refEndNodeIdx + lat_steps); ++endNodeIdx) {
                        
                        if (srcNodeIdx == raceline_index_array[layerIdx] && endNodeIdx == raceline_index_array[dstLayerIdx]) {
                            continue;
                        }

                        Node& endNode = nodesPerLayer[dstLayerIdx][endNodeIdx];
                        
                        MatrixXd path(2, 2);
                        path(0,0) = startNode.x;
                        path(0,1) = startNode.y;
                        path(1,0) = endNode.x;
                        path(1,1) = endNode.y;

                        auto result = calcSplines(path, startNode.psi, endNode.psi);
                        // cout << "result: " << result->el_lengths.size()  << endl;
                        IPair startPoint = make_pair(srcLayerIdx, srcNodeIdx);
                        IPair endPoint = make_pair(dstLayerIdx, endNodeIdx);

                        splineMap[startPoint][endPoint] = *result;

                        // graph에 넣는 과정 
                        wayptGraph.addEdge(startPoint, endPoint);

                        // cout << "startPoint:" << startPoint.first << ", " << startPoint.second << " -> ";
                        // cout << "endPoint:" << endPoint.first << ", " << endPoint.second << endl;
                        
                    }
            }
        }
        return wayptGraph;
    }

    // 단일 스플라인에 대한 spline에 대한 샘플링
    // 샘플링한 점에 대해서 kappa 계산
    auto samplingSpline(MatrixXd &coeffs_x, MatrixXd &coeffs_y, YAML::Node &params) -> pair<vector<Vector2d>, VectorXd> {

        if (coeffs_x.rows() != coeffs_y.rows())
            {
                throw invalid_argument("Coefficient matrices must have the same length!");
            }

        if (coeffs_x.cols() == 2 && coeffs_y.cols() == 2) {
            throw invalid_argument("Coefficient matrices do not have two dimensions!");
        }
        int no_interp_points =  params["sampling"]["no_interp_points"].as<int>();

        VectorXd t_steps(no_interp_points);
        double step = 1.0 / (no_interp_points - 1);
        for (size_t i = 0; i < no_interp_points; ++i) {
            t_steps[i] = i*step;
        }
        
        vector<Vector2d> points_xy;
        VectorXd kappa(no_interp_points+1);
        // kappa.reserve(no_interp_points + 1);

        for (int i = 0; i < no_interp_points; ++i) {
            double t = t_steps(i);
            double t2 = t * t;
            double t3 = t2 * t;

            // 좌표 계산
            double x = coeffs_x(0, 0) + coeffs_x(0, 1) * t + coeffs_x(0, 2) * t2 + coeffs_x(0, 3) * t3;
            double y = coeffs_y(0, 0) + coeffs_y(0, 1) * t + coeffs_y(0, 2) * t2 + coeffs_y(0, 3) * t3;

            // 1차 미분
            double x_d  = coeffs_x(0, 1) + 2 * coeffs_x(0, 2) * t + 3 * coeffs_x(0, 3) * t2;
            double y_d  = coeffs_y(0, 1) + 2 * coeffs_y(0, 2) * t + 3 * coeffs_y(0, 3) * t2;

            // 2차 미분
            double x_dd = 2 * coeffs_x(0, 2) + 6 * coeffs_x(0, 3) * t;
            double y_dd = 2 * coeffs_y(0, 2) + 6 * coeffs_y(0, 3) * t;

            double denom = pow(x_d * x_d + y_d * y_d, 1.5);

            kappa(i) = (x_d * y_dd - y_d * x_dd) / denom;
            points_xy.emplace_back(x, y);
            
        }

        return {points_xy, kappa};
    }

        ///////////////////////////////////////////////////////////////////
        /////////////////////////////제거 과정///////////////////////////////
        ///////////////////////////////////////////////////////////////////

    void pruneEdge(Graph &wayptGraph,
                   NodeMap &nodesPerLayer) {
        int rmv_cnt = 0;

        for (int layerIdx = 0; layerIdx < nodesPerLayer.size(); ++layerIdx) {
            for (int nodeIdx = 0; nodeIdx < nodesPerLayer[layerIdx].size(); ++nodeIdx) {
                IPair srcNodeIdx = make_pair(layerIdx, nodeIdx);

                IPairVector parents;
                bool isParent = wayptGraph.getParentNodes(srcNodeIdx, parents, static_cast<int>(nodesPerLayer.size()));

                if (!isParent) {
                    // cout << layerIdx << ", " << nodeIdx << endl;
                    // cout << "-------" << endl;
                    IPairVector childs;
                    bool isChild = wayptGraph.getChildNodes(srcNodeIdx, childs);
                    if (isChild) {
                        // cout << layerIdx << ", " << nodeIdx << endl;
                        for (auto& child : childs) {
                            // cout << "remove!" << endl;
                            wayptGraph.removeEdge(srcNodeIdx, child, &splineMap, static_cast<int>(nodesPerLayer.size()));
                            rmv_cnt++;
                        }
                    }
                }
            }
        }
        if (rmv_cnt > 0)
            cout << "Removed splines due to isolated nodes: " << rmv_cnt << endl;
    }

    void calcOfflineCost(IVector &raceline_index_array,
                         YAML::Node &params)  {
        if (splineMap.size() <= 0)
        {
            throw invalid_argument("SplineMap's Size is zero!!");
        }

        float w_raceline = params["cost"]["w_raceline"].as<float>();
        float w_raceline_sat = params["cost"]["w_raceline_sat"].as<float>();
        float w_length = params["cost"]["w_length"].as<float>();
        float w_curv_avg = params["cost"]["w_curv_avg"].as<float>();
        float w_curv_peak = params["cost"]["w_curv_peak"].as<float>();
        float lat_resolution = params["lattice"]["lat_resolution"].as<float>();

        for (auto &[startPoint, endPoints] : splineMap)
        {
            for (auto &[endPoint, spline] : endPoints)
            {
                double offline_cost = 0.0;
                int end_layer = endPoint.first;
                int end_node = endPoint.second;

                // 디버깅용
                // cout << "kappa: ";
                // for (int i = 0; i < spline.kappa.size(); ++i) cout << spline.kappa[i] << " ";
                // cout << endl;

                if (end_layer < 0 || end_layer >= raceline_index_array.size())
                {
                    cerr << "[WARNNING] Skipping spline: end_layer=" << end_layer
                         << " out of bounds (0.." << raceline_index_array.size() - 1 << ")\n";
                    continue;
                }

                if (spline.kappa.size() == 0)
                {
                    // cerr << "[WARNNING] Skipping spline: empty curvature data\n";
                    continue;
                }

                double abs_kappa = spline.kappa.array().abs().sum();
                double s_length = spline.el_lengths.sum();
                // cout << "s_length: " << s_length << endl;
                // average curvature
                offline_cost += w_curv_avg * pow(abs_kappa / float(spline.kappa.size()), 2) * s_length;
                // peak curvature
                double max_min = abs(spline.kappa.array().maxCoeff() - spline.kappa.array().minCoeff());
                offline_cost += w_curv_peak * pow(max_min, 2) * s_length;

                // path length
                offline_cost += w_length * s_length;

                // raceline cost

                double raceline_dist = abs(raceline_index_array[end_layer] - end_node) * lat_resolution;
                double raceline_cost = min(w_raceline * s_length * raceline_dist, w_raceline_sat * s_length);

                offline_cost += raceline_cost;

                spline.cost = offline_cost;
                // cout << "(" << startPoint.first << ", " << startPoint.second << ") " << " -> " << "(" << end_layer << ", " << end_node << "): " << offline_cost << endl;
            }
        }
    }
};