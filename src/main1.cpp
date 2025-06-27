#include "graph_planner.hpp"
#include "config.h"

DMap gtpl_map;
DMap sampling_map;

// CSV를 읽어서 DMap으로 변경 
void readDMapFromCSV(const string& pathname, DMap& map) {
    Document csv(pathname, LabelParams(0, -1), SeparatorParams(';'));
    vector<string> labels = csv.GetColumnNames();

    for (const auto& label : labels)
        map[label] = csv.GetColumn<double>(label);
}

// DMap을 CSV에 작성 
void writeDMapToCSV(const string& pathname, DMap& map, char delimiter = ',') {
    ofstream file(pathname);
    if (!file.is_open()) throw runtime_error("Can't open file.");

    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();

    // Header
    size_t i = 0;
    for (const auto& [key, _] : map) {
        file << key;
        if (++i != num_cols) file << delimiter;
    }
    file << '\n';

    // Row map
    for (size_t row = 0; row < num_rows; ++row) {
        size_t j = 0;
        for (const auto& [_, col] : map) {
            file << col[row];
            if (++j != num_cols) file << delimiter;
        }
        file << '\n';
    }

    file.close();
}

// Debug용 함수: map의 columns, rows 개수 print  
void map_size(DMap& map) {
    size_t num_cols = map.size();
    size_t num_rows = map.begin()->second.size();
    cout << "mapsize(" << num_rows << "," << num_cols << ")" << endl;
}

// Dvector를 Map 구조로 추가(연산)
void addDVectorToMap(DMap &map,
                     string attr,
                     const IVector *idx_array = nullptr) {
    size_t len;
    if (idx_array == nullptr) {
        len = map[__x_ref].size();
    } 
    else {
        len = idx_array->size();
    }
    // cout << "attr: "<< attr << " / len:" << len << endl;

    DVector x_out(len), y_out(len);
    string x_label = "x_" + attr;
    string y_label = "y_" + attr;
    
    if (!attr.compare("bound_r")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__width_right][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__width_right][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("bound_l")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] - map[__x_normvec][i] * map[__width_left][i];
            y_out[i] = map[__y_ref][i] - map[__y_normvec][i] * map[__width_left][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    else if (!attr.compare("raceline")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len; ++i) {
            x_out[i] = map[__x_ref][i] + map[__x_normvec][i] * map[__alpha][i];
            y_out[i] = map[__y_ref][i] + map[__y_normvec][i] * map[__alpha][i];
        }

        // x_label = "x_" + attr;
        // y_label = "y_" + attr;
        map[x_label] = x_out;
        map[y_label] = y_out;
    }
    // i번째와 i-1번째 point의 delta_s 계산 
    // delta_s[0] = 0 
    else if (!attr.compare("delta_s")) {
        // cout << "addDVectorToMap:" << attr << endl;
        for (size_t i = 0; i < len - 1; ++i) {
            x_out[i] = map[__s_racetraj][i+1] - map[__s_racetraj][i]; // 마지막 원소는 0
        }
        map[attr] = x_out; 
    }

    // map_size(map);
}

void samplePointsFromRaceline(const DVector& kappa,
                              const DVector& dist,
                              double d_curve,
                              double d_straight,
                              double curve_th,
                              IVector& idx_array) {

    const size_t n = kappa.size();
    double cur_dist = 0.0;
    double next_dist = 0.0;
    double next_dist_min = 0.0;

    for (size_t i = 0; i < n; ++i) {

        // 곡선이면 최소 거리 갱신
        if ((cur_dist + dist[i]) > next_dist_min && fabs(kappa[i]) > curve_th) {
            next_dist = cur_dist;
        }

        // 다음 샘플링 지점 도달
        if ((cur_dist + dist[i]) > next_dist) {
            idx_array.push_back(static_cast<int>(i));

            if (fabs(kappa[i]) < curve_th) {  // 직선 구간
                next_dist += d_straight;
            } else {  // 곡선 구간
                next_dist += d_curve;
            }

            next_dist_min = cur_dist + d_curve;
        }

        cur_dist += dist[i];
    }

    // for (size_t i=0; i < idx_array.size(); ++i) 
    //     cout << idx_array[i] << endl;
    // cout << "size: " << idx_array.size() << endl;
}

double normalizeAngle(double angle) {
    while (angle > M_PI)  angle -= 2.0 * M_PI;
    while (angle < -M_PI) angle += 2.0 * M_PI;
    return angle;
}

void calcHeading(DVector &x_raceline,
                 DVector &y_raceline,
                 DVector &psi) {

    size_t N = x_raceline.size();
    psi.resize(N);

    // 닫힌 회로 가정. 예외 처리 필요
    double dx, dy;
    for (size_t i = 0; i < N; ++i) {
        
        if (i != N -1) {
            dx = x_raceline[i+1] - x_raceline[i];
            dy = y_raceline[i+1] - y_raceline[i];
        } else {
            dx = x_raceline[0] - x_raceline[N - 1];
            dy = y_raceline[0] - y_raceline[N - 1];
        } 
    psi[i] = atan2(dy, dx) - M_PI_2;
        
    normalizeAngle(psi[i]);

    }
    // cout << i<< ": " << psi[i] << endl;
    // cout << psi.size() << endl;

}

void genNode(const DVector& psi_bound_l,
            const DVector& psi_bound_r,
            const double veh_width,
            const float lat_resolution) {

    const size_t N = sampling_map[__alpha].size();
    IVector raceline_index_array;

    // sampling대로, 즉 layer 별로
    for (size_t i = 0; i < N; ++i)    {
        // raceline 기준 인덱스
        int raceline_index = floor(
            (sampling_map[__width_left][i] - veh_width / (2.0 + sampling_map[__alpha][i]) / lat_resolution));
        raceline_index_array.push_back(raceline_index);

        // 알파 보간 시작점
        double s = sampling_map[__alpha][i] - raceline_index * lat_resolution;

        // lateral offset 값 생성
        DVector temp_alphas;
        for (double a = s; a < sampling_map[__width_right][i] - veh_width / 2.0; a += lat_resolution)
            temp_alphas.push_back(a);

        // 위치 계산
        Vector2d ref_pt(sampling_map[__x_ref][i], sampling_map[__y_ref][i]);
        Vector2d norm_vec(sampling_map[__x_normvec][i], sampling_map[__y_normvec][i]);

        vector<Vector2d> pos_list;
        for (double a : temp_alphas)
            pos_list.push_back(ref_pt + norm_vec * a);

        // heading 보간
        DVector psi_interp;

        // 왼쪽
        if (abs(psi_bound_l[i] - sampling_map[__psi][i]) < M_PI) {
            for (int k = 0; k <= raceline_index; ++k)
                psi_interp.push_back(normalizeAngle(psi_bound_l[i] + (sampling_map[__psi][i] - psi_bound_l[i]) * k / raceline_index));
        }
        else {
            double bl = psi_bound_l[i] + 2 * M_PI * (psi_bound_l[i] < 0);
            double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
            for (int k = 0; k <= raceline_index; ++k)
                psi_interp.push_back(normalizeAngle(bl + (p - bl) * k / raceline_index));
        }

        // 오른쪽
        int remain = temp_alphas.size() - raceline_index;
        DVector psi_interp_r;
        if (abs(psi_bound_r[i] - sampling_map[__psi][i]) < M_PI) {
            for (int k = 0; k < remain; ++k)
                psi_interp_r.push_back(normalizeAngle(sampling_map[__psi][i] + (psi_bound_r[i] - sampling_map[__psi][i]) * k / remain));
        }
        else {
            double br = psi_bound_r[i] + 2 * M_PI * (psi_bound_r[i] < 0);
            double p = sampling_map[__psi][i] + 2 * M_PI * (sampling_map[__psi][i] < 0);
            for (int k = 0; k < remain; ++k)
                psi_interp_r.push_back(normalizeAngle(p + (br - p) * k / remain));
        }

        psi_interp.pop_back(); // 중복 제거
        psi_interp.insert(psi_interp.end(), psi_interp_r.begin(), psi_interp_r.end());

        // layer 자료구조 저장

        }
    }

void plotHeading(const DVector &x,
                 const DVector &y,
                 const DVector &psi,
                 double scale = 0.5) {

    double dx, dy;
    double theta, arrow_len;
    double angle;
    double x_arrow1, y_arrow1;
    double x_arrow2, y_arrow2;

    for (size_t i = 0; i < x.size(); ++i) {
        dx = scale * cos(psi[i] + M_PI_2);
        dy = scale * sin(psi[i] + M_PI_2);

        // psi 방향 
        DVector x_line = {x[i], x[i] + dx};
        DVector y_line = {y[i], y[i] + dy};
        plt::plot(x_line, y_line, {{"color", "green"}});

        #if 1
        // 화살촉 
        theta = std::atan2(dy, dx);
        arrow_len = 0.2 * scale;
        angle = M_PI / 6.0;  // 30 degrees

        x_arrow1 = x[i] + dx - arrow_len * cos(theta - angle);
        y_arrow1 = y[i] + dy - arrow_len * sin(theta - angle);

        x_arrow2 = x[i] + dx - arrow_len * cos(theta + angle);
        y_arrow2 = y[i] + dy - arrow_len * sin(theta + angle);

        // 화살촉 그리기 
        plt::plot({x[i] + dx, x_arrow1}, {y[i] + dy, y_arrow1}, {{"color", "green"}});
        plt::plot({x[i] + dx, x_arrow2}, {y[i] + dy, y_arrow2}, {{"color", "green"}});
        #endif

    }
        // raceline 좌표와 psi 프린팅 
        #if 0
        for (size_t i = 0; i < x.size(); ++i) {
        ostringstream label;
        label.precision(2);
        label << fixed << "(" << x[i] << ", " << y[i] << ")\nψ=" << psi[i];

        plt::text(x[i], y[i], label.str());
        #endif
    
}

void visual(const vector<double> &psi_bound_l, const vector<double> &psi_bound_r) {
    plt::clf();

	// plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "black"}});
	// plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "black"}});
    plt::plot(gtpl_map[__x_ref], gtpl_map[__y_ref], {{"color", "blue"}});
    plt::plot(gtpl_map[__x_raceline], gtpl_map[__y_raceline], {{"color", "red"}});
    plt::scatter(sampling_map[__x_raceline], sampling_map[__y_raceline], 30.0, {{"color", "red"}});

    // plotHeading(sampling_map[__x_raceline],
    //             sampling_map[__y_raceline],
    //             sampling_map[__psi]);

    plt::plot(gtpl_map[__x_bound_l], gtpl_map[__y_bound_l], {{"color", "orange"}});
    plt::plot(gtpl_map[__x_bound_r], gtpl_map[__y_bound_r], {{"color", "orange"}});

    // plotHeading(sampling_map[__x_bound_l],
    //             sampling_map[__y_bound_l],
    //             psi_bound_l);

    // plotHeading(sampling_map[__x_bound_r],
    //             sampling_map[__y_bound_r],
    //             psi_bound_r);

    plt::title("Track");
    plt::grid(true);
	plt::axis("equal");
	plt::show();  
}


int main() {
    map<int, int> layer;
    IVector idx_sampling;
    Offline_Params params;

    string map_file_in = "inputs/gtpl_levine.csv";
    string map_file_out = "inputs/gtpl_levine_out.csv";

    // global planner로부터 받은 csv를 기반으로 map에 저장 <label, data> 
    readDMapFromCSV(map_file_in, gtpl_map);

    addDVectorToMap(gtpl_map, "bound_r");
    addDVectorToMap(gtpl_map, "bound_l");
    addDVectorToMap(gtpl_map, "raceline");
    addDVectorToMap(gtpl_map, "delta_s");

    writeDMapToCSV(map_file_out, gtpl_map);

    // layer 간격을 위한 raceline points sampling 
    samplePointsFromRaceline(gtpl_map[__kappa],
                             gtpl_map[__delta_s],
                             params.LON_CURVE_STEP,
                             params.LON_STRAIGHT_STEP,
                             params.CURVE_THR,
                             idx_sampling);

    // cout << "idx size:" << idx_sampling.size() << endl;
    
    for (const auto& [key, vec] : gtpl_map) {
        for (int idx : idx_sampling) {
            if (idx >= 0 && idx < vec.size()) {
                sampling_map[key].push_back(vec[idx]);
            }
        }
    }

    #if 0
    for (int idx : idx_sampling) {
            sampling_map[__x_raceline].push_back(gtpl_map[__x_raceline][idx]);
            sampling_map[__y_raceline].push_back(gtpl_map[__y_raceline][idx]);
            sampling_map[__s_racetraj].push_back(gtpl_map[__s_racetraj][idx]);
            sampling_map[__x_bound_l].push_back(gtpl_map[__x_bound_l][idx]);
            sampling_map[__x_bound_r].push_back(gtpl_map[__x_bound_r][idx]);
            sampling_map[__y_bound_l].push_back(gtpl_map[__y_bound_l][idx]);
            sampling_map[__y_bound_r].push_back(gtpl_map[__y_bound_r][idx]);
        }
    #endif
    // writeDMapToCSV("inputs/sampling_map", sampling_map);
    // map_size(sampling_map); // (51, 3)

    addDVectorToMap(sampling_map, "delta_s", &idx_sampling);
    

    // map_size(sampling_map); // (51, 4)
    writeDMapToCSV("inputs/sampling_mapbefore.csv", sampling_map);
    // 추후 저장될 예정 
    calcHeading(sampling_map[__x_raceline],
                sampling_map[__y_raceline],
                sampling_map[__psi]);
    
    vector<double> psi_bound_l, psi_bound_r;
    
    // 여기서 계산되는 sampling된 bound_l, r은 node 생성 시에만 쓰인다. 
    calcHeading(sampling_map[__x_bound_l],
                sampling_map[__y_bound_l],
                psi_bound_l);

    calcHeading(sampling_map[__x_bound_r],
                sampling_map[__y_bound_r],
                psi_bound_r);   
    
    genNode(psi_bound_l, 
            psi_bound_r,
            params.VEH_WIDTH,
            params.LAT_RESOLUTION);

    // sampling points' info 
    writeDMapToCSV("inputs/sampling_map.csv", sampling_map);
    
    // visual process 
    visual(psi_bound_l, psi_bound_r);
    return 0;
}