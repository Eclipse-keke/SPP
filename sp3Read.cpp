#include"myhead.h"


// ���ļ���ȡ����
std::vector<SatelliteData> readDataFromFile(const std::string& filename, int& timestamp) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "�޷����ļ�: " << filename << std::endl;
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
            // ʱ�����
            timestamp += 5;

        }
        else if (line.rfind("P", 0) == 0) {

            // ����ַ����Ƿ�������ַ��� "PG"
            size_t found = line.find("PG");
            if (found != std::string::npos) {
                // ������
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

// �������ݵ�CSV
void saveDataToCSV(const std::string& filename, int& timestamp, const std::vector<SatelliteData>& data) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "�޷����ļ�: " << filename << std::endl;
        return;
    }

    // д���ͷ
    outfile << "tk,���Ǳ��,X,Y,Z,�Ӳ�" << std::endl;

    // ���������ʽ
    outfile << std::fixed << std::setprecision(6);

    // д������
    for (const SatelliteData& satellite : data) {
        outfile << satellite.timeStamp << "," << satellite.satelliteId << ","
            << satellite.x << "," << satellite.y << ","
            << satellite.z << "," << satellite.clockBias << std::endl;
    }

    outfile.close();
}

