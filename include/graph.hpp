#pragma once
#include "graph_planner.hpp"

class Graph {
private:
    bool isDirected;

public:
    IPairAdjList adjLists;
    explicit Graph(bool directed = true) : isDirected(directed) {}
    
    void addEdge(IPair srcIdx, IPair dstIdx) {
        // adjLists[srcKey].push_back(destIdx);
        adjLists[srcIdx].push_back(dstIdx);
    }

    void printGraph() {
        int size =0;
        for (const auto& [srcNode, childNode] : adjLists) {
            // cout << "(" << srcNode.first << "," << srcNode.second << ")" << ": ";
            for (size_t i = 0; i < childNode.size(); ++i) {
                // cout << "(" << childNode[i].first << ", " << childNode[i].second << ")" << " -> "; 
                size ++;
            }
            // cout << "NULL\n";
            
        }
        // for (const auto& [key, neighbors] : adjLists) {
        //     cout << "(" << get<0>(key) << "," << get<1>(key) << ")" << ": ";
        //     for (int dest : neighbors) {
        //         cout << dest << " -> ";
        //     }
        //     cout << "NULL\n";
        // }
        cout << size << endl;
    }

    bool getChildNodes(const IPair& parentIdx, IPairVector& childIdx) {
        if (adjLists[parentIdx].size() <= 0) 
            return false;
        for (auto& value : adjLists[parentIdx]) {
            childIdx.push_back(value);
            // cout << "("<< value.first << ", " << value.second << ")" << endl;
        }
        return true;
    }

    bool getParentNodes(const IPair& childIdx, IPairVector& parentIdx, int num_layers) {
        parentIdx.clear();
        for (auto &[key, vec] : adjLists) {
            if (childIdx.first == 0) {
                key.first == childIdx.first + num_layers - 1;
                for (const auto &value : vec) {
                    if (value == childIdx) parentIdx.push_back(key);
                    }
            }
            else if (key.first == (childIdx.first) - 1) {
                for (const auto &value : vec) {
                    if (value == childIdx) parentIdx.push_back(key);
                    }
                }
            }
        if (parentIdx.empty()) return false;

        return true;
    }

    void removeEdge(const IPair& srcIdx, const IPair& dstIdx, SplineMap* splineMap, int num_layers) {

        IPairVector& childs = adjLists[srcIdx];
        auto it = remove(childs.begin(), childs.end(), dstIdx);
        if (it != childs.end()) {
            childs.erase(it, childs.end());
        } else {
            return;
        }

        if (splineMap) {
            auto it = splineMap->find(srcIdx);
            if (it != splineMap->end()) {
                auto& splineList = it->second;
                auto it2 = splineList.find(dstIdx);
                if (it2 != splineList.end()) {
                    splineList.erase(it2);
                    if (splineList.empty()) {
                        splineMap->erase(it);
                    }
                }
            }
        }

        if (childs.empty()) {
            IPairVector parents;
            if (getParentNodes(srcIdx, parents, num_layers)) {
                for (const auto& parentIdx : parents) {
                    removeEdge(parentIdx, srcIdx, splineMap, num_layers);
                }
            }
            
        }
    }
};