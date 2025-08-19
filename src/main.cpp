#include "NodeGraph.hpp"
#include "graph_planner.hpp"

DMap loadGlobalTrajectoryMap(string fname);
DMap createSampledTrajectoryMap(DMap gtMap, YAML::Node params);
auto createNodeMap(DMap &stMap, const YAML::Node &params) -> pair<NodeMap, IVector>;
void setInitialPose(DMap &stMap, const NodeMap &nodeMap, const IVector &nodeIndexesOnRaceline, YAML::Node &params);
void visualizeTrajectories(DMap &gtMap, DMap &stMap, const NodeMap &nodeMap, SplineMap &splineMap);

int main() {
    clock_t total_s, total_y;
    total_s = clock();

    string yaml_path = "config/offline_params.yaml";
    YAML::Node params = YAML::LoadFile(yaml_path);

    string map_file_in = "maps/" + params["map_name"].as<string>() + ".csv";
    string map_file_out = "outputs/"+ params["map_name"].as<string>() + "_out.csv";

    DMap gtMap = loadGlobalTrajectoryMap(map_file_in);

    DMap stMap = createSampledTrajectoryMap(gtMap, params);

    auto [nodeMap, nodeIndexesOnRaceline] = createNodeMap(stMap, params);
    // writeDMapToCSV("inputs/stMap.csv", stMap);

    NodeGraph nodeGraph;
    nodeGraph.setNumLayers(nodeMap);
    nodeGraph.genEdges(nodeMap, nodeIndexesOnRaceline, params);
    cout << "Initial generated splines: ";
    nodeGraph.printGraph();

    nodeGraph.writeSplineMapToCSV("outputs/splineMap.csv");

    nodeGraph.pruneEdges(nodeMap, stMap[RL_VX], params);
    cout << "The number of splines generated finally: ";
    nodeGraph.printGraph();

    nodeGraph.computeSplineCost(nodeIndexesOnRaceline, params);
    total_y = clock();

    setInitialPose(stMap, nodeMap, nodeIndexesOnRaceline, params);
    // nodeGraph.printGraph();

    cout << "Total: "<< (double)(total_y - total_s) / CLOCKS_PER_SEC << "s" << endl;

    // printSplineInfo(splineMap, nodeMap);

    visualizeTrajectories(gtMap, stMap, nodeMap, nodeGraph.getSplineMap());

    return 0;
}