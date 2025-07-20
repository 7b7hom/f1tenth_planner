#include "graph_planner.hpp"

Graph::Graph(bool directed) {
    isDirected = directed;
}

void Graph::addEdge(IPair srcNodeIdx, IPair childNodeIdx) {
    // adjLists[srcKey].push_back(destIdx);
    adjLists[srcNodeIdx].push_back(childNodeIdx);
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

void Graph::getChildIdx(IPair srcNodeIdx, IPairVector& childNodeIdx) {
    if (adjLists[srcNodeIdx].size() <= 0) 
        throw runtime_error{"Unable to print child node for srcNodeIdx"};
    for (auto& value : adjLists[srcNodeIdx]) {
        childNodeIdx.push_back(value);
        // cout << "("<< value.first << ", " << value.second << ")" << endl;
    }
}
// 코드 수정 필요. 제기능은 함.
void Graph::getParentNode(IPair &srcNodeIdx, IPairVector &parent) {
    for (auto &[key, vec] : adjLists) {
        if (key.first == (srcNodeIdx.first) - 1) {
            for (const auto &value : vec) {
                if (value == srcNodeIdx) parent.push_back(key);
                }
            }
        }
    }

// adjLists[srcNodeIdx]에서 delNodeIdx만 제거
void Graph::removeEdge(IPair& srcNodeIdx, IPair& delNodeIdx) {
    IPairVector& neighbors = adjLists[srcNodeIdx];
    neighbors.erase(remove(neighbors.begin(), neighbors.end(), delNodeIdx), neighbors.end());
}