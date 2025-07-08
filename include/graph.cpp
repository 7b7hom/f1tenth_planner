#include "graph_planner.hpp"

Graph::Graph(bool directed) {
    isDirected = directed;
}

void Graph::addEdge(IPair srcNodeIdx, IPair childNodeIdx) {
    // adjLists[srcKey].push_back(destIdx);
    adjLists[srcNodeIdx].push_back(childNodeIdx);
}

void Graph::printGraph() {
    for (const auto& [srcNode, childNode] : adjLists) {
        cout << "(" << srcNode.first << "," << srcNode.second << ")" << ": ";
        for (size_t i = 0; i < childNode.size(); ++i) {
            cout << "(" << childNode[i].first << ", " << childNode[i].second << ")" << " -> "; 
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
}

void Graph::getChildIdx(IPair srcNodeIdx, IPairVector& childNodeIdx) {
    if (adjLists[srcNodeIdx].size() <= 0) 
        throw runtime_error{"Unable to print child node for srcNodeIdx"};
    for (const auto& value : adjLists[srcNodeIdx]) {
        childNodeIdx.push_back(value);
    }
}
// 코드 수정 필요. 제기능은 함.  
// void Graph::getParentNode(int target_layer, int value, IPairVector& parent) {
//     for (auto& [key, vec] : adjLists) {
//         if (get<0>(key) == target_layer) {
//                 for (auto it = adjLists[key].begin(); it != adjLists[key].end(); it++) {
//                     if (*it == value) {
//                         parent.push_back(key);
//                     }
//                 }
                    
//             }
//         }
//     }

// void Graph::removeEdge(IPair& parent, int value) {
//     for (auto& [key, vec] : adjLists) {
//         if (key == parent) {
//             auto it = remove(vec.begin(), vec.end(), value);
//             if (it != vec.end()) {
//                 vec.erase(it, vec.end());
//             }
//         }
//     }
// }
