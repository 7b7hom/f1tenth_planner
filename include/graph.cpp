#include "graph_planner.hpp"

Graph::Graph(bool directed) {
    isDirected = directed;
}

void Graph::addEdge(IPair srcIdx, IPair dstIdx) {
    // adjLists[srcKey].push_back(destIdx);
    adjLists[srcIdx].push_back(dstIdx);
}

void Graph::printGraph() {
    int size =0;
    for (const auto& [srcNode, childNode] : adjLists) {
        cout << "(" << srcNode.first << "," << srcNode.second << ")" << ": ";
        for (size_t i = 0; i < childNode.size(); ++i) {
            cout << "(" << childNode[i].first << ", " << childNode[i].second << ")" << " -> "; 
            size ++;
        }
        cout << "NULL\n";
        
    }
    // for (const auto& [key, neighbors] : adjLists) {
    //     cout << "(" << get<0>(key) << "," << get<1>(key) << ")" << ": ";
    //     for (int dest : neighbors) {
    //         cout << dest << " -> ";
    //     }
    //     cout << "NULL\n";
    // }
    cout << "spline 개수" << size << endl;
}

void Graph::getChildNodes(IPair& parentIdx, IPairVector& childIdx) {
    if (adjLists[parentIdx].size() <= 0) 
        throw runtime_error{"Unable to print child node for srcNodeIdx"};
    for (auto& value : adjLists[parentIdx]) {
        childIdx.push_back(value);
        // cout << "("<< value.first << ", " << value.second << ")" << endl;
    }
}
// 코드 수정 필요. 제기능은 함.
void Graph::getParentNodes(IPair& childIdx, IPairVector& parentIdx) {
    for (auto &[key, vec] : adjLists) {
        if (key.first == (childIdx.first) - 1) {
            for (const auto &value : vec) {
                if (value == childIdx) parentIdx.push_back(key);
                }
            }
        }
    }

// adjLists[srcNodeIdx]에서 delNodeIdx만 제거
void Graph::removeEdge(IPair& srcIdx, IPair& dstIdx) {
    IPairVector& neighbors = adjLists[srcIdx];
    neighbors.erase(remove(neighbors.begin(), neighbors.end(), dstIdx), neighbors.end());
}

void Graph::