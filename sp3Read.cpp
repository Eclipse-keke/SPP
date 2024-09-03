#include"myhead.h"


// 从文件读取数据
std::vector<SatelliteData> readDataFromFile(const std::string& filename, int& timestamp) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return {};
    }

    std::string line;
    bool firstStarEnded = false;
    std::vector<SatelliteData> data;

    while (std::getline(infile, line)) {
        if (!firstStarEnded) {
            if (line.find("/* PCV:IGS20      OL/AL:FES2014b NONE     YN ORB:CoN CLK:CoN  ") != std::string::npos) {
                firstStarEnded = true;
            }
            continue;
        }

        if (line.rfind("*", 0) == 0) {
            // 时间戳行
            timestamp += 5;

        }
        else if (line.rfind("P", 0) == 0) {

            // 检查字符串是否包含子字符串 "PG"
            size_t found = line.find("PG");
            if (found != std::string::npos) {
                // 数据行
                SatelliteData satellite;
                satellite.timeStamp = timestamp;
                std::istringstream ss(line);
                ss >> satellite.satelliteId >> satellite.x >> satellite.y >> satellite.z >> satellite.clockBias;
                data.push_back(satellite);
            }
            else {
                continue;
            }
        }
    }

    infile.close();
    return data;
}

// 保存数据到CSV
void saveDataToCSV(const std::string& filename, int& timestamp, const std::vector<SatelliteData>& data) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 写入表头
    outfile << "tk,卫星编号,X,Y,Z,钟差" << std::endl;

    // 设置输出格式
    outfile << std::fixed << std::setprecision(6);

    // 写入数据
    for (const SatelliteData& satellite : data) {
        outfile << satellite.timeStamp << "," << satellite.satelliteId << ","
            << satellite.x << "," << satellite.y << ","
            << satellite.z << "," << satellite.clockBias << std::endl;
    }

    outfile.close();
}

