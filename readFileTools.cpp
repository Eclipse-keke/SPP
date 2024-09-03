#define _CRT_SECURE_NO_WARNINGS

#include"myhead.h"
void setstr(char* des, const char* src, int n)
{
	char* p = des;
	const char* q = src;
	while (*q && q < src + n)
	{
		*p++ = *q++;
	}
	*p-- = '\0';
	//去掉尾部空格
	while (p >= des && *p == ' ')
	{
		*p-- = '\0';
	}
}
/*
获取文件数据的行数，是从END OF HEADER之后开始起算的
返回值：row是最后要返回的行数，当读到END OF HEADER时，flag = 1，row = 1，row开始计数。
*/
int getrow(FILE* fp_nav)
{
	char line[100];// Assuming maximum line length is 100 characters
	int row = 0;
	int flag = 0;

	//Seek to the begining of the file
	fseek(fp_nav, 0, SEEK_SET);

	//Read the file line by line
	//function fgets can read the content of fp_nav and export to line
	while (fgets(line, sizeof(line), fp_nav) != NULL)
	{
		//Check if the line contatins "END of Header"
		if (strstr(line, "END OF HEADER") != NULL)
		{
			flag = 1;
			row = 1;  //Start counting rows from the line after "END OF HEADER"
			continue;
		}

		//If "END OF HEADER"has been found,start counting rows
		if (flag == 1)
		{
			row++;
		}
	}
	return row;

}


/*
这个函数
buff，读进来的东西存在buff里
i,从第列开始读
n,读多少个
*/
static double strtonum(const std::string& buff, int i, int n) {
	double value = 0.0;
	char str[256] = { 0 };
	char* p = str;

	/*****************
	*  当出现以下三种情况报错，返回0.0
	* 1.起始位置<0
	* 2.读取字符串个数<i
	* 3.str里面存放的字节数<n
	****************/

	if (i < 0 || buff.length() < static_cast<size_t>(i) || sizeof(str) - 1 < n) {
		return 0.0;
	}

	for (size_t j = i; j < buff.length() && n > 0; ++j, --n) {
		//三目操作符：D和d为文件中科学计数法部分，将其转换成二进制能读懂的e
		*p++ = (buff[j] == 'D' || buff[j] == 'd') ? 'e' : buff[j];
	}
	*p = '\0';
	//三目操作符，将str中存放的数以格式化读取到value中。
	return sscanf(str, "%lf", &value) == 1 ? value : 0.0;
}

//将字符串转换为浮点数,i起始位置，n输入多少个字符
static double strtonum(const char* buff, int i, int n)
{
	double value = 0.0;
	char str[256] = { 0 };
	char* p = str;
	/************************************
	* 当出现以下三种情况报错，返回0.0
	* 1.起始位置<0
	* 2.读取字符串个数<i
	* 3.str里面存放的字节数<n
	*************************************/
	if (i < 0 || (int)strlen(buff) < i || (int)sizeof(str) - 1 < n)
	{
		return 0.0;
	}
	for (buff += i; *buff && --n >= 0; buff++)
	{
		//三目操作符：D和d为文件中科学计数法部分，将其转换成二进制能读懂的e
		*p++ = ((*buff == 'D' || *buff == 'd') ? 'e' : *buff);
	}
	*p = '\0';
	//三目操作符，将str中存放的数以格式化读取到value中。
	return sscanf(str, "%lf", &value) == 1 ? value : 0.0;
}
/*
读取N文件
输入参数：
fp_nav：导航文件头文件指针
nav_h:头文件结构体
*/
void readrinex_head(FILE* fp_nav, pnav_head nav_h)
{
	char buff[MAXRINEX] = { 0 };
	char* label = buff;
	int i = 0;
	int j = 0;
	while (fgets(buff, MAXRINEX, fp_nav))
	{
		if (strstr(label, "RINEX VERSION / TYPE"))
		{
			nav_h->ver = strtonum(buff, 0, 9);
			strncpy((nav_h->type), buff + 20, 16);
			continue;
		} 
		else if (strstr(label, "GPSA"))
		{
			nav_h->ION_alpha[0] = strtonum(buff, 6, 12);
			nav_h->ION_alpha[1] = strtonum(buff, 6 + 12, 12);
			nav_h->ION_alpha[2] = strtonum(buff, 6 + 12 + 12, 12);
			nav_h->ION_alpha[3] = strtonum(buff, 6 + 12 + 12 + 12, 12);
			continue;
		}
		else if (strstr(label, "GPSB"))
		{
			nav_h->ION_beta[0] = strtonum(buff, 6, 12);
			nav_h->ION_beta[1] = strtonum(buff, 6 + 12, 12);
			nav_h->ION_beta[2] = strtonum(buff, 6 + 12 + 12, 12);
			nav_h->ION_beta[3] = strtonum(buff, 6 + 12 + 12 + 12, 12);
			continue;
		}
	}
	std::cout << "这里已经读完了头文件" << std::endl;

}


/*
dataBlockRead
This function is for Reading Data Block
if choose is C, read Beidou
if choose is G, read GPS
*/
void dataBlocksRead(std::vector<std::vector<std::string>>dataBlocks, std::vector<nav_body>&  nav_b_vector, char choose)
{

	for (std::vector<std::string>& datablock : dataBlocks)
	{
		nav_body nav_b;//这里已经初始化了
		//Read the First Line
		//std::cout << datablock[0] << std::endl;
		nav_b.monument = datablock[0].at(0);
		nav_b.sPRN = (int)strtonum(datablock[0], 1, 2);
		nav_b.TOC_Y = (int)strtonum(datablock[0], 4, 4);
		nav_b.TOC_M = (int)strtonum(datablock[0], 9, 2);
		nav_b.TOC_D = (int)strtonum(datablock[0], 12, 2);
		nav_b.TOC_H = (int)strtonum(datablock[0], 15, 2);
		nav_b.TOC_Min = (int)strtonum(datablock[0], 18, 2);
		nav_b.TOC_Sec = (int)strtonum(datablock[0], 21, 2);
		nav_b.sa0 = strtonum(datablock[0], 23, 19);
		nav_b.sa1 = strtonum(datablock[0], 23 + 19, 19);
		nav_b.sa2 = strtonum(datablock[0], 23 + 19 + 19, 19);
		//Read the Second Line
		nav_b.IODE = strtonum(datablock[1], 4, 19);
		nav_b.Crs = strtonum(datablock[1], 4 + 19, 19);
		nav_b.delta_n = strtonum(datablock[1], 4 + 19 + 19, 19);
		nav_b.M0 = strtonum(datablock[1], 4 + 19 + 19 + 19, 19);

		// Read the Third Line
		nav_b.Cuc = strtonum(datablock[2], 4, 19);
		nav_b.e = strtonum(datablock[2], 4 + 19, 19);
		nav_b.Cus = strtonum(datablock[2], 4 + 19 + 19, 19);
		nav_b.sqrtA = strtonum(datablock[2], 4 + 19 + 19 + 19, 19);

		// Read the Fourth Line
		nav_b.TOE = strtonum(datablock[3], 4, 19);
		nav_b.Cic = strtonum(datablock[3], 4 + 19, 19);
		nav_b.OMEGA = strtonum(datablock[3], 4 + 19 + 19, 19);
		nav_b.Cis = strtonum(datablock[3], 4 + 19 + 19 + 19, 19);

		// Read the Fifth Line
		nav_b.i0 = strtonum(datablock[4], 4, 19);
		nav_b.Crc = strtonum(datablock[4], 4 + 19, 19);
		nav_b.omega = strtonum(datablock[4], 4 + 19 + 19, 19);
		nav_b.deltaomega = strtonum(datablock[4], 4 + 19 + 19 + 19, 19);

		// Read the Sixth Line
		nav_b.IDOT = strtonum(datablock[5], 4, 19);
		nav_b.L2code = strtonum(datablock[5], 4 + 19, 19);
		nav_b.GPSweek = strtonum(datablock[5], 4 + 19 + 19, 19);
		nav_b.L2Pflag = strtonum(datablock[5], 4 + 19 + 19 + 19, 19);

		// Read the 7th Line
		nav_b.sACC = strtonum(datablock[6], 4, 19);
		nav_b.sHEA = strtonum(datablock[6], 4 + 19, 19);
		nav_b.TGD = strtonum(datablock[6], 4 + 19 + 19, 19);
		nav_b.IODC = strtonum(datablock[6], 4 + 19 + 19 + 19, 19);

		// Read the Last line
		nav_b.TTN = strtonum(datablock[7], 4, 19);
		nav_b.fit = strtonum(datablock[7], 4 + 19, 19);
		nav_b.spare1 = strtonum(datablock[7], 4 + 19 + 19, 19);
		nav_b.spare2 = strtonum(datablock[7], 4 + 19 + 19 + 19, 19);

		if (nav_b.monument == choose)
			nav_b_vector.push_back(nav_b);
	}
}

/*
This function is for reading the body of the rinex file
input:
const char* fileName
n_n: The number of lines except for the head file.
return:
0: successfully.
1:Error opening file.
*/
int readrinex_body(const char* fileName, std::vector<nav_body>& nav_b_vector, int n_n,char choose)
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
		}
	}

	if (!foundEndOfHeader) {
		std::cerr << "Error: 'END OF HEADER' not found!" << std::endl;
		return 1;
	}
	
	std::vector<std::vector<std::string>> dataBlocks;//It is for storing a lot of data blocks.
	for (int i = 0; i < n_n / 8; i++)//n_n/8就是数据块的数量（一个数据块有8行）
	{
		std::vector<std::string> dataBlock;
		int linesAfterEndHeader = 0;
		// 输出"END OF HEADER"后的内容
		while (std::getline(file, line)) {
			//std::cout << line << std::endl;
			dataBlock.push_back(line);
			linesAfterEndHeader++;
			if (linesAfterEndHeader >= 8)
			{
				dataBlocks.push_back(dataBlock);
				break;
			}
		}
	}

	dataBlocksRead(dataBlocks, nav_b_vector,choose);
	file.close();
	return 0;
}

