#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "include/qcustomplot.h"
#include <QApplication>
#include <QMainWindow>

using namespace std;

bool readXYcsv(const string& filename, vector<double>& xs, vector<double>& ys)
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "파일을 열 수 없습니다: " << filename << endl;
        return false;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string x_str, y_str;
        if (!getline(ss, x_str, ',')) continue;
        if (!getline(ss, y_str, ',')) continue;

        xs.push_back(stod(x_str));
        ys.push_back(stod(y_str));
    }
    return true;
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    vector<double> x_left, y_left;
    vector<double> x_right, y_right;

    if (!readXYcsv("bound_l.csv", x_left, y_left)) return -1;
    if (!readXYcsv("bound_r.csv", x_right, y_right)) return -1;

    QMainWindow window;
    QCustomPlot *customPlot = new QCustomPlot();

    // 왼쪽 경계선 그래프 추가
    customPlot->addGraph();
    customPlot->graph(0)->setPen(QPen(Qt::blue));  // 파란색
    customPlot->graph(0)->setData(QVector<double>::fromStdVector(x_left),
                                 QVector<double>::fromStdVector(y_left));

    // 오른쪽 경계선 그래프 추가
    customPlot->addGraph();
    customPlot->graph(1)->setPen(QPen(Qt::red));  // 빨간색
    customPlot->graph(1)->setData(QVector<double>::fromStdVector(x_right),
                                 QVector<double>::fromStdVector(y_right));

    customPlot->xAxis->setLabel("X [m]");
    customPlot->yAxis->setLabel("Y [m]");
    customPlot->rescaleAxes();
    customPlot->replot();

    window.setCentralWidget(customPlot);
    window.resize(800, 600);
    window.show();

    return app.exec();
}
