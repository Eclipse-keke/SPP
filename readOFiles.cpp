#include"myhead.h"


// �Զ��庯�������ڴ��ַ����л�ȡָ�����ȵ����ַ�������ʹ�����ǿյ�Ҳ�ܶ���ָ�������Ŀո�
std::string get_fixed_length_substr(const std::string& str, size_t start, size_t length) {
	std::string result;
	size_t actualLength = length < str.size() - start ? length: (str.size() - start);

	// ����ָ�����ȵ��ַ�
	result = str.substr(start, actualLength);

	// ���ʵ�ʳ���С�����賤�ȣ�����ո�
	if (actualLength < length) 
	{
		result.append(length - actualLength, ' ');
	}

	return result;
}

OData readOFile(const std::string& r_o_name)
{
	std::ifstream file(r_o_name);
	if (!file.is_open())
	{
		throw std::runtime_error("�޷����ļ���");
	}

	OData rinexData;
	//"ABMF00GLP_R_20240250000_01D_30S_MO.rnx"�����ֿ�ʼ��ȡ
	rinexData.station = r_o_name.substr(0, r_o_name.size() - 4);
	std::cout << "��ȡ���Ĳ�վ��Ϊ��" << rinexData.station << std::endl;

	int epoch_cnt = -1;
	std::vector<double> epoch_buf;//�����Ԫ
	std::vector<std::vector<double>> data_buf;
	std::vector<int> index_buf;//�������
	std::vector<std::vector<double>> time_buf;
	int read_state = -1;

	std::string line;
	int counter = 0;
	while (std::getline(file, line))
	{
		if (line.find("TIME OF FIRST OBS") != std::string::npos)
		{
			//  2024     1    25     0     0    0.0000000     GPS         TIME OF FIRST OBS����ȡ2024     1    25     0     0    0.0000000�Ĳ���
			std::string dateStr = line.substr(0, line.size() - 35);
			std::istringstream dateStream(dateStr);
			//�Ѷ�ȡ����Ϣ����rinexData
			dateStream >> rinexData.date.year >> rinexData.date.month >> rinexData.date.day >> rinexData.date.hour >> rinexData.date.minute >> rinexData.date.second;
		}

		if (line.find("APPROX POSITION XYZ") != std::string::npos)
		{
			std::string rcvposStr = line.substr(0, line.size() - 22);
			//��������  2919785.7120 -5383745.0670  1774604.6920
			//���Ľ��Ʊ��λ��(��λ:�ף��Ƽ�ʹ��ITRSϵͳ)�ƶ�ƽ̨��ѡ
			std::istringstream rcvposStream(rcvposStr);
			rcvposCoor value;
			while (rcvposStream >> value.x >> value.y >> value.z)
			{
				rinexData.rcvpos.push_back(value);//istringStream���Զ������м�����Ŀո�
			}
		}
		if (line.find("SYS / # / OBS TYPES") != std::string::npos)
		{
			char sysChar = line[0];
			if (sysChar == 'G')//�����GPS���ǣ��Ǿ������У���Ҫ�����¶�һ��
			{
				std::string line2;
				if (std::getline(file, line2) && line2.find("SYS / # / OBS TYPES") != std::string::npos)
				{
					//G   18 C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L
					//           S2L C5Q L5Q D5Q S5Q                                  SYS / # / OBS TYPES
					int numTypes = std::stoi(line.substr(4, 2));
					std::string typeStr = line.substr(7, 13 * 4);//����13��������5��
					if (!line2.empty())
					{
						typeStr += line2.substr(0, line2.find("SYS / # / OBS TYPES"));
					}
					std::istringstream iss(typeStr);
					std::string token;
					while (iss >> token)
					{
						rinexData.type.push_back(token);
					}
				}
			}
		}
		//��ʼ��ʱ����> 2024 01 25 00 00  0.0000000  0 47
//year (4 digits), month, day, hour, min(2 digits)
//sec
//Epoch flag: 0��ʾOK, 
//1��ʾpower failure between previous and current epoch
//>1��ʾ�������
//47��ʾ�۲⵽��47������
		if (line[0] == '>')
		{
			read_state = 1; //��ʶ�Ѿ�������ʱ���ʶ�У�������ȡ״̬
			std::vector<double> timeValues;
			std::istringstream timeStream(line.substr(1));
			double value;
			while (timeStream >> value)
			{
				timeValues.push_back(value);
			}
			time_buf.push_back(timeValues);

			//������ΪO�ļ�30���¼һ�Σ�����epoch_cntֻ����0��30��60Ȼ��һֱ���ϼ�
			if (epoch_cnt == -1)
			{
				epoch_cnt = 0;
			}
			else
			{
				epoch_cnt += 30;
			}
		}

		if (read_state == 1 && line[0] != '>')
		{
			char sysChar = line[0];
			//����ֻ��GPS���ǵ����ݣ�Ҳ����G����
			std::string tmpStr = get_fixed_length_substr(line, 5, 100);
			std::string segment_C1C = tmpStr.substr(0, 13);
			std::string segment_C2X = tmpStr.substr(64, 13);
			//[]�ķ����ű�ʾ�������κα�����
			//char c��ʾ����һ������c
			//all_off�㷨�Ὣ�ַ�����ÿһ���ַ����ݸ����lambda���ʽ
			//
			bool allSpaces_C1C = std::all_of(segment_C1C.begin(), segment_C1C.end(), [](char c)->bool
			{
				return std::isspace(static_cast<unsigned char>(c));
			});
			bool allSpaces_C2X = std::all_of(segment_C2X.begin(), segment_C2X.end(), [](char c)->bool
				{
					return std::isspace(static_cast<unsigned char>(c));
				});
			if (sysChar == 'G' && !allSpaces_C1C && !allSpaces_C2X)
			{
				//��¼��ʱ��ʱ�䣨�룩
				epoch_buf.push_back(epoch_cnt);
				//��¼���Ǳ��
				index_buf.push_back(std::stoi(line.substr(1, 2)));
				//��ʼ��ȡ����һ��ѵ�����
				std::vector<double> tmp(2, std::nan(" "));//��nanֵ��ʼ��һ��vector.
				//��¼C1C����
				tmp[0] = std::stod(segment_C1C);
				//��¼C2X����
				tmp[1] = std::stod(segment_C2X);

				data_buf.push_back(tmp);
			}
		}
		counter++;
		//std::cout << "�Ѿ���¼��" << counter << "��" << std::endl;
	}
	rinexData.epoch = epoch_buf;
	rinexData.index = index_buf;
	rinexData.data = data_buf;
	rinexData.time = time_buf;
	file.close();
	return rinexData;
}

int readOFileContent(const char* fileName)
{
	std::ifstream file(fileName);

	if (!file.is_open())
	{
		std::cerr << "Error opening file!" << std::endl;
	}

	std::string line;
	bool foundEndOfHeader = false;

	//Read the content of file line by line.
	//getline means read a line from STD input stream and stored it to ''line''.
	while (std::getline(file, line))
	{
		if (line.find("END OF HEADER") != std::string::npos)
		{
			foundEndOfHeader = true;
			std::cout << "�ҵ���END OF HEADER!" << std::endl;
			break;
			return 0;
		}
	}

	if (!foundEndOfHeader) 
	{
		std::cerr << "Error: 'END OF HEADER' not found!" << std::endl;
		return 1;
	}
	std::cout << "Content after 'END OF HEADER':" << std::endl;
	return 0;
}
