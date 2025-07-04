#include "graph_planner.hpp"

Graph::Graph(bool directed) {
    isDirected = directed;
}

void Graph::addEdge(ITuple srcKey, int destIdx) {
    adjLists[srcKey].push_back(destIdx);
}

void Graph::printGraph() {
    for (const auto& [key, neighbors] : adjLists) {
        cout << "(" << get<0>(key) << "," << get<1>(key) << ")" << ": ";
        for (int dest : neighbors) {
            cout << dest << " -> ";
        }
        cout << "NULL\n";
    }
}

void Graph::getChildIdx(ITuple srcKey, IVector& childIdx) {
    if (adjLists[srcKey].size() <= 0) 
        throw runtime_error{"Unable to print child node for srcKey"};
    for (const auto& value : adjLists[srcKey]) {
        childIdx.push_back(value);
    }
}
// 코드 수정 필요. 제기능은 함.  
void Graph::getParentNode(int target_layer, int value, vector<ITuple>& parent) {
    for (auto& [key, vec] : adjLists) {
        if (get<0>(key) == target_layer) {
                for (auto it = adjLists[key].begin(); it != adjLists[key].end(); it++) {
                    if (*it == value) {
                        parent.push_back(key);
                    }
                }
                    
            }
        }
    }

void Graph::removeEdge(ITuple& parent, int value) {
    for (auto& [key, vec] : adjLists) {
        if (key == parent) {
            auto it = remove(vec.begin(), vec.end(), value);
            if (it != vec.end()) {
                vec.erase(it, vec.end());
            }
        }
    }
}
