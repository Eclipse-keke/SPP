#include"myhead.h"


// 自定义函数，用于从字符串中获取指定长度的子字符串，即使后面是空的也能读入指定数量的空格
std::string get_fixed_length_substr(const std::string& str, size_t start, size_t length) {
	std::string result;
	size_t actualLength = length < str.size() - start ? length: (str.size() - start);

	// 复制指定长度的字符
	result = str.substr(start, actualLength);

	// 如果实际长度小于所需长度，补充空格
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
		throw std::runtime_error("无法打开文件！");
	}

	OData rinexData;
	//"ABMF00GLP_R_20240250000_01D_30S_MO.rnx"对名字开始读取
	rinexData.station = r_o_name.substr(0, r_o_name.size() - 4);
	std::cout << "读取到的测站名为：" << rinexData.station << std::endl;

	int epoch_cnt = -1;
	std::vector<double> epoch_buf;//存放历元
	std::vector<std::vector<double>> data_buf;
	std::vector<int> index_buf;//存放索引
	std::vector<std::vector<double>> time_buf;
	int read_state = -1;

	std::string line;
	int counter = 0;
	while (std::getline(file, line))
	{
		if (line.find("TIME OF FIRST OBS") != std::string::npos)
		{
			//  2024     1    25     0     0    0.0000000     GPS         TIME OF FIRST OBS，读取2024     1    25     0     0    0.0000000的部分
			std::string dateStr = line.substr(0, line.size() - 35);
			std::istringstream dateStream(dateStr);
			//把读取的信息存入rinexData
			dateStream >> rinexData.date.year >> rinexData.date.month >> rinexData.date.day >> rinexData.date.hour >> rinexData.date.minute >> rinexData.date.second;
		}

		if (line.find("APPROX POSITION XYZ") != std::string::npos)
		{
			std::string rcvposStr = line.substr(0, line.size() - 22);
			//读入数据  2919785.7120 -5383745.0670  1774604.6920
			//地心近似标记位置(单位:米，推荐使用ITRS系统)移动平台可选
			std::istringstream rcvposStream(rcvposStr);
			rcvposCoor value;
			while (rcvposStream >> value.x >> value.y >> value.z)
			{
				rinexData.rcvpos.push_back(value);//istringStream能自动忽略中间隔开的空格
			}
		}
		if (line.find("SYS / # / OBS TYPES") != std::string::npos)
		{
			char sysChar = line[0];
			if (sysChar == 'G')//如果是GPS卫星，那就有两行，需要再往下读一行
			{
				std::string line2;
				if (std::getline(file, line2) && line2.find("SYS / # / OBS TYPES") != std::string::npos)
				{
					//G   18 C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L
					//           S2L C5Q L5Q D5Q S5Q                                  SYS / # / OBS TYPES
					int numTypes = std::stoi(line.substr(4, 2));
					std::string typeStr = line.substr(7, 13 * 4);//读了13个，还有5个
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
		//开始读时间行> 2024 01 25 00 00  0.0000000  0 47
//year (4 digits), month, day, hour, min(2 digits)
//sec
//Epoch flag: 0表示OK, 
//1表示power failure between previous and current epoch
//>1表示特殊情况
//47表示观测到了47颗卫星
		if (line[0] == '>')
		{
			read_state = 1; //标识已经读到了时间标识行，开启读取状态
			std::vector<double> timeValues;
			std::istringstream timeStream(line.substr(1));
			double value;
			while (timeStream >> value)
			{
				timeValues.push_back(value);
			}
			time_buf.push_back(timeValues);

			//这是因为O文件30秒记录一次，所以epoch_cnt只能是0、30、60然后一直往上加
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
			//我们只读GPS卫星的数据，也就是G数据
			std::string tmpStr = get_fixed_length_substr(line, 5, 100);
			std::string segment_C1C = tmpStr.substr(0, 13);
			std::string segment_C2X = tmpStr.substr(64, 13);
			//[]的方括号表示不捕获任何变量。
			//char c表示接受一个参数c
			//all_off算法会将字符串中每一个字符传递给这个lambda表达式
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
				//记录此时的时间（秒）
				epoch_buf.push_back(epoch_cnt);
				//记录卫星编号
				index_buf.push_back(std::stoi(line.substr(1, 2)));
				//开始读取下面一大堆的数据
				std::vector<double> tmp(2, std::nan(" "));//用nan值初始化一个vector.
				//记录C1C数据
				tmp[0] = std::stod(segment_C1C);
				//记录C2X数据
				tmp[1] = std::stod(segment_C2X);

				data_buf.push_back(tmp);
			}
		}
		counter++;
		//std::cout << "已经记录了" << counter << "行" << std::endl;
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
			std::cout << "找到了END OF HEADER!" << std::endl;
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
