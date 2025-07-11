#include "graph_planner.hpp"
#include "config.h"

// gen_spline.cpp에 정의된 함수 선언
void runPlanningPipeline(const Offline_Params& params, const std::string& map_file_in, const std::string& map_file_out);

int main() {
    Offline_Params params;
    std::string map_file_in = "/home/uiiiqns/f1tenth_planner/inputs/gtpl_levine.csv";
    std::string map_file_out = "/home/uiiiqns/f1tenth_planner/inputs/gtpl_levine_out.csv";

    runPlanningPipeline(params, map_file_in, map_file_out);

    return 0;
}