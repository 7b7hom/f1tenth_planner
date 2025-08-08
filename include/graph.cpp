#include "graph_planner.hpp"

Graph::Graph(bool directed) {
    isDirected = directed;
}

// 아직 spline 계수나 길이 확정 X 단순히 두 노드 사이에 edge가 존재한다는 사실만을 기록할 때
void Graph::addEdge(ITuple srcKey, int destIdx) {
    adjLists[srcKey][destIdx] = EdgeInfo(); 
}

// edge에 대한 모든 상세 정보가 확정되었을 때 호출
void Graph::addEdge(ITuple srcKey, int destIdx, const Eigen::RowVector4d& coeffs_x, const Eigen::RowVector4d& coeffs_y, double spline_len) {
    EdgeInfo info;
    info.spline_len = spline_len;
    info.coeffs_x_orig = coeffs_x;
    info.coeffs_y_orig = coeffs_y;
    info.offline_cost = 0.0; // 초기 비용은 0.0

    adjLists[srcKey][destIdx] = info;
}

void Graph::printGraph() {
    for (const auto& [srcKey, destMap] : adjLists) {
        std::cout << "(" << std::get<0>(srcKey) << "," << std::get<1>(srcKey) << ")" << ": ";

        double total_outgoing_cost_for_node = 0.0;
        if (!destMap.empty()) {
            for (const auto& [destIdx, edgeInfo] : destMap) {
                total_outgoing_cost_for_node += edgeInfo.offline_cost;
            }
        }

        if (destMap.empty()) {
            std::cout << "NULL";
        } else {
            for (const auto& [destIdx, _] : destMap) { 
                std::cout << destIdx << " -> ";
            }
            std::cout << "NULL";
        }

        std::cout << " (" << std::fixed << std::setprecision(2) << total_outgoing_cost_for_node << ")\n"; 
    }
}

void Graph::getChildIdx(ITuple srcKey, IVector& childIdx) {
    childIdx.clear();
    auto it = adjLists.find(srcKey);
    if (it != adjLists.end()) {
        for (const auto& [destIdx, _] : it->second) {
            childIdx.push_back(destIdx);
        }
    } else {
        throw std::runtime_error("Node with srcKey not found in adjLists for getChildIdx.");
    }
}

void Graph::getParentNode(int target_layer, int value, std::vector<ITuple>& parent, int num_layers) {
    parent.clear();

    int expected_parent_layer;
    if (target_layer == 0) {
        expected_parent_layer = num_layers - 1; 
    } else {
        expected_parent_layer = target_layer - 1; 
    }

    for (auto& [key, destMap] : adjLists) {
        if (std::get<0>(key) == expected_parent_layer) {
            if (destMap.count(value)) { 
                parent.push_back(key);
            }
        }
    }
}

void Graph::removeEdge(const ITuple& parent, int value) {
    auto it_parent_map = adjLists.find(parent);
    if (it_parent_map != adjLists.end()) {
        it_parent_map->second.erase(value);
    }
}

// 읽기 전용
const GraphAdjListsMap& Graph::getAdjLists() const {
    return adjLists;
}

// 쓰기 전용
GraphAdjListsMap& Graph::getAdjLists_mutable() {
    return adjLists;
}
