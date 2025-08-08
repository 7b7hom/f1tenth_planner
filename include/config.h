#pragma once

#include <string>

const int rl_kappa = 9;
const int rl_s = 7;

struct Offline_Params {

    // [lattice]
    float lat_resolution = 0.5;         // lateral spacing (m) between nodes along each normal
    float variable_heading = true;      // determines if node headings are interpolated between boundary and race line angles (else: match race line)
    
    float lon_straight_step = 4.0;      // max distance (m) between norm vectors along straights on the reference line
    float lon_curve_step = 0.9;         // max norm vector spacing (m) on reference-line curves

    float curve_thr = 0.21;             // recommend: 0.15 ~ 0.3
    float lat_offset = 2.0;

    int max_lat_steps = 3;              // permitted lateral deviation from the raceline per meter traveled
    float virt_goal_n = true;           // proxy target node per layer
    float min_vel_race = 0.0;           // minimum allowed speed as % of global race line
    float max_lateral_accel = 10.0;
    float closure_detection_dist = 0.5; // if track ends(first, last) are within this distance (m), treat as closed loop

    // [planningtarget]
    float vel_decrease_lat = 0.1;       // planning target velocity: % goal speed reduced per meter the goal is offset from the raceline
    float min_plan_horizon = 300;       // minimum number of layers (or distance) included in the online planning graph
    std::string plan_horizon_mode = "distance"; // set planning horizon mode

    // [sampling]
    float stepsize_approx = 2.5;        // spline sampling interval

    // [vehicle]
    float veh_width = 1.0;              // vehicle width (m)
    float veh_length = 4.7;             // vehicle length in m
    float veh_turn = 5.0;               // min turn radius (m)

    // [cost]
    float w_raceline = 1.0;             // penalty for path length and lateral offset from raceline
    float w_raceline_sat = 1.0;         // max race line cost per meter due to lateral offset
    float w_length = 0.0;               // penalty applied for each meter of spline path
    float w_curv_avg = 7500.0;          // penalty factor for average curvature
    float w_curv_peak = 2500.0;         // penalty factor for highest curvature
    float w_virt_goal = 10000.0;        // penalty per meter of lateral offset at the virtual goal node

    // [custom]
    float max_heading_offset = M_PI / 4;
};

