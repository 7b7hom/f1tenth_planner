#pragma once
#include "graph_planner.hpp"

class NodeGraph {
private:
  IPairAdjList nodeGraph;
  SplineMap splineMap;
  int num_layers;

public:
  void setNumLayers(NodeMap& nodeMap) {
      num_layers = static_cast<int>(nodeMap.size());
  }

  SplineMap &getSplineMap() { return splineMap; }

  Spline &at(const IPair &start, const IPair &end)  {
    return splineMap[start][end];
  }

  const Spline &at(const IPair &start, const IPair &end) const  {
    return splineMap.at(start).at(end);
  }

  void addEdge(IPair srcIdx, IPair dstIdx)  {
    nodeGraph[srcIdx].push_back(dstIdx);
  }

  void printGraph()  {
    int size = 0;
    for (const auto &[srcNode, childNode] : nodeGraph)    {
      // cout << "(" << srcNode.first << "," << srcNode.second << ")" << ": ";
      for (size_t i = 0; i < childNode.size(); ++i)      {
        // cout << "(" << childNode[i].first << ", " << childNode[i].second << ")" << " -> ";
        size++;
      }
      // cout << "NULL\n";
    }
    // for (const auto& [key, neighbors] : nodeGraph) {
    //     cout << "(" << get<0>(key) << "," << get<1>(key) << ")" << ": ";
    //     for (int dest : neighbors) {
    //         cout << dest << " -> ";
    //     }
    //     cout << "NULL\n";
    // }
    cout << size << endl;
  }

  IPairVector getChildList(const IPair &srcIdx)  {
    IPairVector childList;
    for (auto &value : nodeGraph[srcIdx])    {
      childList.push_back(value);
      // cout << "("<< value.first << ", " << value.second << ")" << endl;
    }
    return childList;
  }

  IPairVector getParentList(const IPair &srcIdx)  {
    IPairVector parentList;
    for (auto &[key, vec] : nodeGraph)    {
      if (srcIdx.first == 0)      {
        key.first == srcIdx.first + num_layers - 1;
        for (const auto &value : vec)        {
          if (value == srcIdx)
            parentList.push_back(key);
        }
      }
      else if (key.first == (srcIdx.first) - 1)      {
        for (const auto &value : vec)
        {
          if (value == srcIdx)
            parentList.push_back(key);
        }
      }
    }
    return parentList; // 길이가 0인경우 밖에서 처리.
  }

  void removeEdge(const IPair &srcIdx, const IPair &dstIdx)  {
      IPairVector &childs = nodeGraph[srcIdx];
      auto it = remove(childs.begin(), childs.end(), dstIdx);
      if (it != childs.end()) {
          childs.erase(it, childs.end());
      } else {
          return;
      }

      auto itSpline = splineMap.find(srcIdx);
      if (itSpline != splineMap.end()) {
          auto &splineList = itSpline->second;
          auto it2 = splineList.find(dstIdx);
          if (it2 != splineList.end()) {
              splineList.erase(it2);
              if (splineList.empty()) {
                  splineMap.erase(itSpline);
              }
          }
      }

      if (childs.empty()) {
          IPairVector parentList = getParentList(srcIdx);
          if (!parentList.empty()) {
              for (const auto &parentIdx : parentList) {
                  removeEdge(parentIdx, srcIdx);
              }
          }
      }
  }

  void writeSplineMapToCSV(const string &filename)  {
    ofstream fout(filename);
    if (!fout.is_open())
      throw runtime_error("Cannot open file");

    // 헤더
    fout << "start_layer,start_idx,end_layer,end_idx,coeffs_x,coeffs_y,kappa,points_x,points_y,raceline\n";

    for (const auto &[start, endMap] : splineMap)    {
      for (const auto &[end, spline] : endMap)      {
        fout << start.first << "," << start.second << ","
             << end.first << "," << end.second << ",";

        // coeffs_x
        for (int i = 0; i < spline.coeffs_x.rows(); ++i)
          for (int j = 0; j < spline.coeffs_x.cols(); ++j)
            fout << spline.coeffs_x(i, j) << (i == spline.coeffs_x.rows() - 1 && j == spline.coeffs_x.cols() - 1 ? "," : " ");

        // coeffs_y
        for (int i = 0; i < spline.coeffs_y.rows(); ++i)
          for (int j = 0; j < spline.coeffs_y.cols(); ++j)
            fout << spline.coeffs_y(i, j) << (i == spline.coeffs_y.rows() - 1 && j == spline.coeffs_y.cols() - 1 ? "," : " ");

        // kappa
        for (int i = 0; i < spline.kappaVector.size(); ++i)
          fout << spline.kappaVector(i) << (i == spline.kappaVector.size() - 1 ? "," : " ");

        // points_xy
        for (int i = 0; i < spline.points_xy.size(); ++i)
          fout << spline.points_xy[i].x() << " " << spline.points_xy[i].y()
               << (i == spline.points_xy.size() - 1 ? "," : " ");

        fout << (spline.raceline ? 1 : 0) << "\n";
      }
    }
    fout.close();
  }

  void readSplineMapFromCSV(const string &filename)
  {
    ifstream fin(filename);
    if (!fin.is_open())
      throw runtime_error("Cannot open file");

    string line;
    getline(fin, line); // 헤더 스킵

    while (getline(fin, line))    {
      stringstream ss(line);
      string item;

      IPair start, end;
      Spline spline;

      // start_layer, start_idx, end_layer, end_idx
      getline(ss, item, ',');
      start.first = stoi(item);
      getline(ss, item, ',');
      start.second = stoi(item);
      getline(ss, item, ',');
      end.first = stoi(item);
      getline(ss, item, ',');
      end.second = stoi(item);

      // coeffs_x
      getline(ss, item, ',');
      stringstream sx(item);
      DVector vx;
      double val;
      while (sx >> val)
        vx.push_back(val);

      spline.coeffs_x = MatrixXd(4, 1);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 1; j++)
          spline.coeffs_x(i, j) = vx[i * 1 + j];

      // coeffs_y
      getline(ss, item, ',');
      stringstream sy(item);
      DVector vy;
      while (sy >> val)
        vy.push_back(val);
      spline.coeffs_y = MatrixXd(4, 1);
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 1; j++)
          spline.coeffs_y(i, j) = vy[i * 1 + j];

      // kappa
      getline(ss, item, ',');
      stringstream sk(item);
      DVector vk;
      while (sk >> val)
        vk.push_back(val);
      spline.kappaVector = VectorXd::Map(vk.data(), vk.size());

      // points_xy
      getline(ss, item, ',');
      stringstream sp(item);
      vector<Vector2d> pts;
      double x, y;
      while (sp >> x >> y)
        pts.emplace_back(x, y);
      spline.points_xy = pts;

      // raceline
      getline(ss, item, ',');
      spline.raceline = (stoi(item) != 0);

      splineMap[start][end] = spline;
    }

    fin.close();
  }

  void genEdges(NodeMap &nodeMap,
                const IVector &nodeIndexesOnRaceline,
                YAML::Node &params) {

    float lat_offset = params["lattice"]["lat_offset"].as<float>();
    float lat_resolution = params["lattice"]["lat_resolution"].as<float>();
    int max_lat_steps = params["lattice"]["max_lat_steps"].as<int>();

    if (lat_offset <= 0.0)    {
      throw invalid_argument("Too small lateral offset!");
    }

    // raceline
    for (int layerIdx = 0; layerIdx < num_layers; ++layerIdx)    {
      int dstLayerIdx = layerIdx + 1;
      if (dstLayerIdx >= num_layers)      {
        dstLayerIdx -= num_layers;
      }
      int startNodeIdx = nodeIndexesOnRaceline[layerIdx];
      int endNodeIdx = nodeIndexesOnRaceline[dstLayerIdx];

      Node &startNode = nodeMap[layerIdx][startNodeIdx];
      Node &endNode = nodeMap[dstLayerIdx][endNodeIdx];

      MatrixXd path(2, 2);
      path(0, 0) = startNode.x;
      path(0, 1) = startNode.y;
      path(1, 0) = endNode.x;
      path(1, 1) = endNode.y;

      auto result = computeSplines(path, startNode.psi, endNode.psi);

      IPair startPoint = make_pair(layerIdx, startNodeIdx);
      IPair endPoint = make_pair(dstLayerIdx, endNodeIdx);

      result->raceline = true;

      splineMap[startPoint][endPoint] = *result;

      addEdge(startPoint, endPoint);
    }

    // cout << num_layers << endl; // 출력: 51
    // 레이어 별 loop
    // raceline spline 먼저 생성해서 splineMap에 등록, 뒤에서 재등록하지 않게 index 겹치면 pass
    for (int layerIdx = 0; layerIdx < num_layers; ++layerIdx)    {

      int srcLayerIdx = layerIdx;
      int dstLayerIdx = layerIdx + 1;

      // cout << "srcLayerIdx:" << srcLayerIdx << endl;
      // cout << "num_layers" << num_layers << endl;

      // 마지막 layer의 경우 0번째 layer와 연결시킬 수 있도록 dstLayerIdx 조정
      if (dstLayerIdx >= num_layers)      {
        dstLayerIdx -= num_layers;
      }

      // start layer 내 노드별 loop
      for (size_t srcNodeIdx = 0; srcNodeIdx < nodeMap[srcLayerIdx].size(); ++srcNodeIdx)      {
        // 기준 노드
        Node &startNode = nodeMap[srcLayerIdx][srcNodeIdx];

        int refEndNodeIdx = nodeIndexesOnRaceline[dstLayerIdx] - (nodeIndexesOnRaceline[srcLayerIdx] - srcNodeIdx);
        refEndNodeIdx = max(0, min(refEndNodeIdx, static_cast<int>(nodeMap[dstLayerIdx].size() - 1)));
        // refEndNodeIdx = clamp(refEndNodeIdx, 0, static_cast<int>(nodeMap[dstLayerIdx].size() - 1));
        // int refEndNodeIdx = srcNodeIdx;
        // int refEndNodeIdx = nodeIndexesOnRaceline[dstLayerIdx];
        Node &endNode = nodeMap[dstLayerIdx][refEndNodeIdx];

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
             endNodeIdx <= min(static_cast<int>(nodeMap[dstLayerIdx].size() - 1), refEndNodeIdx + lat_steps); ++endNodeIdx)        {

          if (srcNodeIdx == nodeIndexesOnRaceline[layerIdx] && endNodeIdx == nodeIndexesOnRaceline[dstLayerIdx])          {
            continue;
          }

          Node &endNode = nodeMap[dstLayerIdx][endNodeIdx];

          MatrixXd path(2, 2);
          path(0, 0) = startNode.x;
          path(0, 1) = startNode.y;
          path(1, 0) = endNode.x;
          path(1, 1) = endNode.y;

          auto result = computeSplines(path, startNode.psi, endNode.psi);
          // cout << "result: " << result->el_lengths.size()  << endl;
          IPair startPoint = make_pair(srcLayerIdx, srcNodeIdx);
          IPair endPoint = make_pair(dstLayerIdx, endNodeIdx);

          splineMap[startPoint][endPoint] = *result;

          // graph에 넣는 과정
          addEdge(startPoint, endPoint);

          // cout << "startPoint:" << startPoint.first << ", " << startPoint.second << " -> ";
          // cout << "endPoint:" << endPoint.first << ", " << endPoint.second << endl;
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////
  /////////////////////////////제거 과정///////////////////////////////
  ///////////////////////////////////////////////////////////////////

  void pruneEdges(NodeMap &nodeMap, const DVector& raceline_vx, YAML::Node &params)  {

    float veh_turn = params["vehicle"]["veh_turn"].as<float>();
    float min_vel_race = params["lattice"]["min_vel_race"].as<float>();
    float max_lateral_accel = params["lattice"]["max_lateral_accel"].as<float>();
    float veh_width = params["vehicle"]["veh_width"].as<float>();
    int rmv_cnt = 0;
    
    // s_time = clock();

    for (size_t layer_idx = 0; layer_idx < num_layers;++layer_idx) {
      int srcLayerIdx = layer_idx;
      for (size_t node_idx = 0; node_idx < nodeMap[srcLayerIdx].size(); ++node_idx) {

        IPair start = make_pair(layer_idx ,node_idx);
        IPairVector childList = getChildList(start);

        // 연결된 노드와 loop
        for (auto& end : childList) {
            MatrixXd& coeffs_x = splineMap[start][end].coeffs_x;
            MatrixXd& coeffs_y = splineMap[start][end].coeffs_y;
            // spline 위의 점들을 샘플링(no_interp_points개수만큼)
            auto [points_xy, kappaVector] = sampleSingleSpline(coeffs_x, coeffs_y, params);
            // 점들을 기준으로 pruneEdges에 가서 1.곡률 2.트랙내 여부 에 따라 remove를 한다.
            if (kappaVector.size() == 0 || points_xy.size() == 0) {
                cerr << "Invalid spline sampling" << endl;
                continue;
            }
            // 해당 spline위에서 샘플링한 점이 track을 벗어나면 pruneEdges()에 갈 수 있도록.
            splineMap[start][end].kappaVector = kappaVector;
            splineMap[start][end].points_xy = points_xy;

            if (!splineMap[start][end].raceline) {

                int layer_idx = start.first;
                double vel_rl = raceline_vx[layer_idx] * min_vel_race;
                double min_turn = pow(vel_rl, 2) / max_lateral_accel;

                bool toRemove = false;

                for (int j = 0; j < kappaVector.size(); ++j) {

                    double kappa_val = abs(kappaVector(j));
                    
                    if ((kappa_val > 1.0 / veh_turn || kappa_val > 1.0 / min_turn))
                    {
                        toRemove = true;
                        break;
                    }
                }
                if (toRemove) {
                    removeEdge(start, end);
                    rmv_cnt++;
                    }
                }
            }
        }
    }

    for (int layerIdx = 0; layerIdx < num_layers; ++layerIdx) {
      for (int nodeIdx = 0; nodeIdx < nodeMap[layerIdx].size(); ++nodeIdx) {
        IPair srcNodeIdx = make_pair(layerIdx, nodeIdx);

        IPairVector parentList = getParentList(srcNodeIdx);

        if (parentList.empty())        {
          // cout << layerIdx << ", " << nodeIdx << endl;
          // cout << "-------" << endl;
          IPairVector childList = getChildList(srcNodeIdx);
          if (!childList.empty())          {
            // cout << layerIdx << ", " << nodeIdx << endl;
            for (auto &child : childList)            {
              // cout << "remove!" << endl;
              removeEdge(srcNodeIdx, child);
              rmv_cnt++;
            }
          }
        }
      }
    }
    if (rmv_cnt > 0)
      cout << "Removed splines due to curvature conditions && isolated nodes: " << rmv_cnt << endl;
  }

  void computeSplineCost(IVector &nodeIndexesOnRaceline,
                       YAML::Node &params)  {
    if (splineMap.size() <= 0)    {
      throw invalid_argument("SplineMap's Size is zero!!");
    }

    float w_raceline = params["cost"]["w_raceline"].as<float>();
    float w_raceline_sat = params["cost"]["w_raceline_sat"].as<float>();
    float w_length = params["cost"]["w_length"].as<float>();
    float w_curv_avg = params["cost"]["w_curv_avg"].as<float>();
    float w_curv_peak = params["cost"]["w_curv_peak"].as<float>();
    float lat_resolution = params["lattice"]["lat_resolution"].as<float>();

    for (auto &[startPoint, endPoints] : splineMap)    {
      for (auto &[endPoint, spline] : endPoints)      {
        double offline_cost = 0.0;
        int end_layer = endPoint.first;
        int end_node = endPoint.second;

        // 디버깅용
        // cout << "kappa: ";
        // for (int i = 0; i < spline.kappa.size(); ++i) cout << spline.kappa[i] << " ";
        // cout << endl;

        if (end_layer < 0 || end_layer >= nodeIndexesOnRaceline.size())        {
          cerr << "[WARNNING] Skipping spline: end_layer=" << end_layer
               << " out of bounds (0.." << nodeIndexesOnRaceline.size() - 1 << ")\n";
          continue;
        }

        if (spline.kappaVector.size() == 0)        {          
          // cerr << "[WARNNING] Skipping spline: empty curvature data\n";
          continue;
        }

        double abs_kappa = spline.kappaVector.array().abs().sum();
        double s_length = spline.el_lengths.sum();
        // cout << "s_length: " << s_length << endl;

        // average curvature
        offline_cost += w_curv_avg * pow(abs_kappa / float(spline.kappaVector.size()), 2) * s_length;
        // peak curvature
        double max_min = abs(spline.kappaVector.array().maxCoeff() - spline.kappaVector.array().minCoeff());
        offline_cost += w_curv_peak * pow(max_min, 2) * s_length;

        // path length
        offline_cost += w_length * s_length;

        // raceline cost
        double raceline_dist = abs(nodeIndexesOnRaceline[end_layer] - end_node) * lat_resolution;
        double raceline_cost = min(w_raceline * s_length * raceline_dist, w_raceline_sat * s_length);

        offline_cost += raceline_cost;

        spline.cost = offline_cost;
        // cout << "(" << startPoint.first << ", " << startPoint.second << ") " << " -> " << "(" << end_layer << ", " << end_node << "): " << offline_cost << endl;
      }
    }
  }
};