#pragma once
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>  
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "rapidcsv.h"
#include "matplotlibcpp.h"

#define __x_ref "x_ref_m"
#define __y_ref " y_ref_m"
#define __width_right " width_right_m"
#define __width_left " width_left_m"
#define __x_normvec " x_normvec_m"
#define __y_normvec " y_normvec_m"
#define __alpha " alpha_m"
#define __kappa " kappa_racetraj_radpm"
#define __s_racetraj " s_racetraj_m"
#define __psi " psi_racetraj_rad"

#define __x_raceline "x_raceline"
#define __y_raceline "y_raceline"
#define __x_bound_r "x_bound_r"
#define __y_bound_r "y_bound_r"
#define __x_bound_l "x_bound_l"
#define __y_bound_l "y_bound_l"
#define __delta_s "delta_s"
#define __psi_bound_l "psi_bound_l"
#define __psi_bound_r "psi_bound_r"

using namespace std;
using namespace rapidcsv;
using namespace Eigen;
namespace plt = matplotlibcpp;

// 스플라인 결과를 담기 위한 구조체: x, y 방향 계수, 행렬 M, 정규화된 노멀 벡터
struct SplineResult {
    MatrixXd coeffs_x;              // 각 구간의 x 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd coeffs_y;              // 각 구간의 y 방향 3차 다항식 계수 행렬 (구간 개수 x 4)
    MatrixXd M;                     // 스플라인 계수 계산에 사용된 시스템 행렬
    MatrixXd normvec_normalized;    // 각 구간의 법선 벡터를 정규화한 값 (구간 개수 x 2)
};

struct Node {
    int layer_idx;
    int node_idx;
    double x;
    double y;
    double psi;
    double kappa;
    bool raceline;
};

typedef vector<double> DVector;
typedef vector<int>    IVector;
typedef map<string, DVector> DMap;
typedef map<string, IVector> IMap;
typedef vector<vector<Node>> NodeMap;
typedef tuple<int, int> ITuple;
typedef map<ITuple, IVector> TupleMap;

// TUM의 GraphBase 역할 
class Graph {
private:
    TupleMap adjList;
    bool isDirected;
public:
    Graph(bool directed = true) {
        isDirected = directed;
    }
    
    void addEdge(ITuple srcKey, int destIdx) {
        adjList[srcKey].push_back(destIdx);
    }

    void printGraph() const {
        for (const auto& [key, neighbors] : adjList) {
            cout << "(" << get<0>(key) << "," << get<1>(key) << ")" << ": ";
            for (int dest : neighbors) {
                cout << dest << " -> ";
            }
            cout << "NULL\n";
        }
    }

    void getChildIdx(ITuple srcKey, IVector& childIdx) {
        if (adjList[srcKey].size() <= 0) 
            throw runtime_error{"Unable to print child node for srcKey"};
        for (const auto& value : adjList[srcKey]) {
            childIdx.push_back(value);
        }
    }
    // 코드 수정 필요. 제기능은 함.  
    void getParentNode(int target_layer, int value, vector<ITuple>& parent) {
        for (auto& [key, vec] : adjList) {
            if (get<0>(key) == target_layer) {
                    for (auto it = adjList[key].begin(); it != adjList[key].end(); it++) {
                        if (*it == value) {
                            parent.push_back(key);
                        }
                    }
                        
                }
            }
        }

    void removeEdge(ITuple& parent, int value) {
        for (auto& [key, vec] : adjList) {
            if (key == parent) {
                auto it = remove(vec.begin(), vec.end(), value);
                if (it != vec.end()) {
                    vec.erase(it, vec.end());
                }
            }
        }
    }

};


