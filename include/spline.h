#pragma once
#include "graph_planner.hpp"
#include "graph.h"

unique_ptr<Spline> calcSplines(const MatrixXd &path,
                               double psi_s,
                               double psi_e,
                               bool use_dist_scaling = true);
VectorXd calcKappa(MatrixXd &coeffs_x,
                   MatrixXd &coeffs_y,
                   VectorXd &t_steps);
pair<Graph, SplineMap> genEdges(DMap &stMap,
                                NodeMap &nodesPerLayer,
                                const IVector &raceline_index_array,
                                YAML::Node &params);

pair<VectorXd, VectorXd> interpSplines(DMap &stMap,
                                       MatrixXd &coeffs_x,
                                       MatrixXd &coeffs_y,
                                       float stepsize_approx,
                                       const float &veh_width,
                                       double spline_len = NAN,
                                       int no_interp_points = 10);
bool checkInsideBounds(DMap &stMap,
                       const Vector2d &pos,
                       const float veh_width);