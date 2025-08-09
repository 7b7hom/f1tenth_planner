#pragma once
#include "graph_planner.hpp"

class Graph {
private:
    bool isDirected;
public:
    IPairAdjList adjLists;
    Graph(bool directed = true);
    void addEdge(IPair srcIdx, IPair dstIdx);
    void printGraph();
    bool getChildNodes(const IPair& parentIdx, IPairVector& childIdx);
    bool getParentNodes(const IPair& childIdx, IPairVector& parentIdx, int num_layers);
    void removeEdge(const IPair& srcIdx, const IPair& dstIdx, SplineMap* splineMap, int& remove_cnt, int num_layers);
};