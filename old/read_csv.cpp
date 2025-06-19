#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ifstream file;
    file.open("tgtpl_levine.csv", ios_base::in);
    if(!file.is_open()) {
        cout << "Error: 파일이 열리지 않음." << endl;
        return 1;
    }

    string line;
    while (getline(file, line)) {
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        
        if (line.empty() || line[0] == '#') continue;
        
            stringstream ss(line); // line을 문자 stream으로 변환 
            string value;
            vector<double> numbers;

            while (getline(ss, value, ';')) {
                try {
                    numbers.push_back(stod(value));
                } catch (...) {
                    std::cout << "숫자 변환 실패: [" << value << "]" << std::endl;
                }
                    
            }
        cout << "x: " << numbers[0] << ", y: " << numbers[1] <<", width_right: " << numbers[2]<< endl;
    
    }
    file.close();
    return 0;
}