#pragma once

#include <string>

const int rl_kappa = 9;
const int rl_s = 7;

struct Offline_Params {

        // [LATTICE]
        float LAT_RESOLUTION=0.16;       // Lateral spacing (m) between nodes along each normal
        float VARIABLE_HEADING=true;    // Determines if node headings are interpolated between boundary and race line angles (else: match race line).
        
        float LON_STRAIGHT_STEP=2.0;   // Max distance (m) between norm vectors along straights on the reference line
        float LON_CURVE_STEP=0.9;      // Max norm vector spacing (m) on reference-line curves.
                                        
        float CURVE_THR=0.25;          // recommand: 0.15 ~ 0.3
        float LAT_OFFSET=0.45;

        int MAX_LAT_STEPS=3;    // Permitted lateral deviation from the raceline per meter traveled.
        float VIRT_GOAL_N=true;         // Proxy target node per layer. Graph search uses this node directly, skipping iterative checks on adjacent goal-layer nodes.
        float MIN_VEL_RACE=0.5;        // Minimum allowed speed as % of global race line; splines below this are excluded (set to 0.0 to allow all).
        float MAX_LATERAL_ACCEL=4.0;
        float CLOSURE_DETECTION_DIST=0.1; // If track ends(first, last) are within this distance (m), treat as closed loop.
        
        // [PLANNINGTARGET]
        float VEL_DECREASE_LAT=0.1;     // PLANNING TARGET VELOCITY : % goal speed reduced per meter the goal is offset from the raceline.
        float MIN_PLAN_HORIZON=300;     // Minimum number of layers (or distance) included in the online planning graph.
        std::string PLAN_HORIZON_MODE="distance";   // Set planning horizon mode: distance-based or fixed layers
                                                    // Note: For 'layers', adjust lon_straight_step and lon_curve_step.

        // [SAMPLING]
        float STEPSIZE_APPROX=2.5;  // Spline sampling interval (smaller values improve collision-check accuracy)

        // [VEHICLE]
        float VEH_WIDTH=0.4;    // Vehicle width (m); increase for greater safety margin
        float VEH_LENGTH=4.7;   // vehicle length in m (in order to obtain a safety distance in follow mode)
        float VEH_TURN=3.0;     // Min turn radius (m); splines below this are excluded.

        // [COST]
        float W_RACELINE=1.0;   // Penalty for path length and lateral offset from raceline, applied per meter
        float W_RACELINE_SAT=1.0;   // Max race line cost per meter due to lateral offset.
        float W_LENGTH=0.0;     // Penalty applied for each meter of spline path.
        float W_CURV_AVG=7500.0;    // Penalty factor for average curvature along each meter of path.
        float W_CURV_PEAK=2500.0;   // Penalty factor for highest curvature along the spline, applied per meter.
        float W_VIRT_GOAL=10000.0;  // Penalty per meter of lateral offset from the race line at the virtual goal node

};

 
// struct Online_Params {

// };

// class Offline_Config {
// public:
//     const Offline_Params& get() const { return params_; }
//     void set(const Offline_Params& new_params) { params_ = new_params; }

// private:
//     Offline_Params params_;
// };

// class Online_Config {
// public:
//     const Online_Params& get() const { return params_; }
//     void set(const Online_Params& new_params) { params_ = new_params; }

// private:
//     Online_Params params_;
// };

