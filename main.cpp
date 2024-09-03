#define _CRT_SECURE_NO_WARNINGS
#include"myhead.h"
using namespace std;


void saveToCSV(const std::vector<NavData>& navData,const std::string& filename)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "无法打开文件: " << filename << std::endl;
		return;
	}

	// 写入CSV表头
	file << "卫星标识,PRN,tk,X,Y,Z,clockBias" << std::endl;

	for (const auto& nav : navData) {
		file << nav.monument << ","
			<< nav.sPRN << ","
			<< nav.tk<<","
			<< std::setprecision(9) << nav.X << ","
			<< std::setprecision(9) << nav.Y << ","
			<< std::setprecision(9) << nav.Z <<","
			<<std::setprecision(9)<<nav.clockBias<<std::endl;

	}

	file.close();

}



/*
产生t的机器
*/
vector<double> tVectorGenerate()
{
	vector<double> tVector;
	double d240219 = 1;//因为是星期一所以是1
	for (int i = 0; i <= 1440; i += 5)
	{
		double t = d240219 * 86400 + i*60;
		tVector.push_back(t);
	}
	return tVector;
}


/*
功能：找到t-toc最小的那一个星历
输入t（是tVector产生的），输入卫星编号，输入那个存放这个卫星数据的容器
输出那个对应的星历（卫星编号，时间戳）
*/
nav_body findTrueOne(double t,vector<nav_body> vessel)
{
	double delta_t = 9999999;
	double delta_t_min = 9999999;
	int min_index = 0;
	vector<double> delta_ts;
	vector<int>nav_b_indexes;
	for (int i = 0; i < vessel.size(); i++)
	{
		nav_body& nav_b = vessel[i];
		//toc = toc_compute(nav_b.TOE, nav_b.TOC_H, nav_b.TOC_Min, nav_b.TOC_Sec);
		delta_t =fabs( t - nav_b.TOE);
		delta_ts.push_back(delta_t);
		nav_b_indexes.push_back(i); // 记录 nav_b 的索引
		if (delta_t < delta_t_min) // 找到更小的值
		{
			delta_t_min = delta_t;
			min_index = i;
		}
	}

	return vessel[min_index];

}
 


pnav_head readNavDataHead(const char* fileName)
{
	FILE* fp_nav = NULL;
	pnav_head nav_h = NULL;
	fp_nav = fopen(fileName, "r");
	nav_h = (pnav_head)malloc(sizeof(nav_head));//给N文件头开辟空间
	readrinex_head(fp_nav, nav_h);
	cout << "The 8 ionospheric parameters obtained are:" << endl;
	int i = 0;
	for (i = 0; i < 4; i++)
	{
		cout << nav_h->ION_alpha[i] << endl;
	}

	for (i = 0; i < 4; i++)
	{
		cout << nav_h->ION_beta[i] << endl;
	}
	return nav_h;
}




//文件读取函数：
//参数：文件绝对 / 相对路径
//读取N文件的GPS卫星数据，读得的结果是很多个nav body然后存储在一个vector里
vector<nav_body> NDataRead_G(const char* fileName)
{
	//Read Head File
	FILE* fp_nav = NULL;//导航星历文件指针
	pnav_body nav_b = NULL;

	//N文件开始读取
	fp_nav = fopen(fileName, "r");
	int n = getrow(fp_nav); //获取导航星历文件的行数
	std::cout << "导航星历文件的行数是" << endl;
	cout << n << endl;
	rewind(fp_nav);//Move the file pointer to the starting position

	nav_b = (pnav_body)malloc(sizeof(nav_body));//malloc space for body of File N.

	fclose(fp_nav);

	//GPS卫星总数据
	vector<nav_body> g_vector;
	//Read Body File.
	readrinex_body(fileName, g_vector, n,'G');
	return g_vector;
}

//参数：文件绝对 / 相对路径
//读取N文件的Beidou卫星数据，读得的结果是很多个nav body然后存储在一个c_vector里
vector<nav_body> NDataRead_C(const char* fileName)
{
	//Read Head File
	FILE* fp_nav = NULL;//导航星历文件指针
	pnav_body nav_b = NULL;

	//N文件开始读取
	fp_nav = fopen(fileName, "r");
	int n = getrow(fp_nav); //获取导航星历文件的行数
	std::cout << "导航星历文件的行数是" << endl;
	cout << n << endl;
	rewind(fp_nav);//Move the file pointer to the starting position

	nav_b = (pnav_body)malloc(sizeof(nav_body));//malloc space for body of File N.

	fclose(fp_nav);

	//Beidou卫星总数据
	vector<nav_body> c_vector;
	readrinex_body(fileName, c_vector, n,'C');
	return c_vector;
}



/*
工作函数
*/
vector<NavData> myReadRNX(vector<nav_body>g_vector,vector<double>tVector)
{
	//GPS卫星单颗卫星数据，分32个存储
	vector<vector<nav_body>>  gpsSats;
	gpsSats.resize(32);


	cout << "完毕" << endl;
	//将时间转化为gps时
	double tk = 0;
	std::vector<NavData> data;

	for (nav_body& nav_b : g_vector)
	{
		gpsSats[nav_b.sPRN-1].push_back(nav_b);
	}

	for (double& t : tVector)
	{
		for (int i = 0; i < 32; i++)
		{
			nav_body trueParams;
			//1. 对tVector里的每一个t，寻找里这个t最近的toe对应的一组参数，也就是，t-toe>而且t-toe最小
			trueParams = findTrueOne(t, gpsSats[i]);
			//2.用这个t计算，计算t-toe
			NavData navData;
			tk = t - trueParams.TOE;
			navData = myComputeLocationOfSatellite(trueParams, tk);//t是周内秒，一周以内的从星期一86400*7
			navData.time = t-86400;
			data.push_back(navData);
		}
	}
	saveToCSV(data, "brdcResult.csv");
	std::cout << "数据已保存到 brdcResult.csv" << std::endl;
	return data;
}

/*
工作函数2：For Beidou
*/
vector<NavData> myReadRNX_C(vector<nav_body>c_vector, vector<double>tVector)
{
	//北斗卫星单颗卫星数据，分62个存储
	vector<vector<nav_body>>  bdsSats;
	bdsSats.resize(62);


	cout << "完毕" << endl;
	//将时间转化为gps时
	double tk = 0;
	std::vector<NavData> data;

	for (nav_body& nav_b : c_vector)
	{
		bdsSats[nav_b.sPRN - 1].push_back(nav_b);
	}

	for (double& t : tVector)
	{
		for (int i = 0; i < 62; i++)
		{

			if (bdsSats[i].size() == 0)
				continue;
			nav_body trueParams;
			//1. 对tVector里的每一个t，寻找里这个t最近的toe对应的一组参数，也就是，t-toe>而且t-toe最小
			trueParams = findTrueOne(t, bdsSats[i]);
			//2.用这个t计算，计算t-toe
			NavData navData;
			tk = t - trueParams.TOE-14;//北斗卫星要减去14
			navData = myComputeLocationOfSatellite(trueParams, tk);//t是周内秒，一周以内的从星期一86400*7
			navData.time = t - 86400;
			data.push_back(navData);
		}
	}
	saveToCSV(data, "brdcResult_C.csv");
	std::cout << "数据已保存到 brdcResult_C.csv" << std::endl;
	return data;
}


/*
卫星精密星历计算函数
输入：卫星参数结构体，时间t
输出：卫星在瞬时地固坐标系下的坐标
*/
NavData myComputeLocationOfSatellite(nav_body nav_b,double tk)
{
	double n = 0;
	double M = 0;
	double E = 0;
	double f = 0;
	double u_pie = 0;
	double delta_u = 0;
	double delta_r = 0;
	double delta_i = 0;
	double u = 0;
	double r = 0;
	double i = 0;
	double x = 0;
	double y = 0;
	double L = 0;
	double X = 0;
	double Y = 0;
	double Z = 0;
	double clockBias = 0;
	n = computeAverageAngularVelocity(nav_b.sqrtA,  nav_b.delta_n);
	M = computeInstantMeanAnomaly(nav_b.M0,  n, tk);
	E = computeEccentricAnomaly(M, nav_b.e);
	f = solveTrueAnomaly(E, nav_b.e);
	u_pie = solveArgumentOfLatitude(nav_b.omega, f);
	solvePerturbationCorrection(nav_b.Cuc, nav_b.Cus, nav_b.Crc, nav_b.Crs, nav_b.Cic, nav_b.Cis, &delta_u, &delta_r, &delta_i, u_pie);
	perturbationCorrection(nav_b.sqrtA, &u, &r, &i, nav_b.IDOT, E, tk, nav_b.e, u_pie, nav_b.i0, delta_u, delta_r, delta_i);
	computexyInOrbitPlane(r, u, &x, &y);
	L = computeMomentOMEGA(tk, nav_b.deltaomega, nav_b.OMEGA, nav_b.TOE);
	solveSatillitePosition(L, i, x, y, &X, &Y, &Z);
	clockBias = clockBiasCompute(nav_b, tk);
	/*
	cout << "年月日时分秒"<<nav_b.TOC_Y<<" "<<nav_b.TOC_M<<" "<<nav_b.TOC_D<<" "<<nav_b.TOC_H<<" "<<nav_b.TOC_Min<<" "<<nav_b.TOC_Sec << endl;
	cout << "卫星标识:" << nav_b.monument << endl;
	cout << "PRN是" << nav_b.sPRN << endl;
	cout << "计算得到的卫星坐标是" << endl;
	cout << "X:"<<std::setprecision(9) << X << endl;
	cout << "Y:" << std::setprecision(9) << Y << endl;
	cout << "Z:" << std::setprecision(9) << Z << endl;
	cout << "钟差是" << endl;
	cout << "clockBias:" << std::setprecision(9) << clockBias <<endl;
	*/
	NavData navData;
	navData.monument = nav_b.monument;
	navData.sPRN = nav_b.sPRN;
	navData.tk = tk;
	navData.X = X;
	navData.Y = Y;
	navData.Z = Z;
	navData.clockBias = clockBias;
	return navData;
}

//用来读取sp3的函数
// 输入文件名 inputFilename
int myReadSp3(std::string inputFilename) {
	std::string outputFilename = "sp3output.csv"; // 输出文件名
	int timestamp = -5;
	std::vector<SatelliteData> data = readDataFromFile(inputFilename, timestamp);
	if (!data.empty()) {
		saveDataToCSV(outputFilename, timestamp, data);
		std::cout << "数据已保存到 " << outputFilename << std::endl;
	}
	else {
		std::cerr << "没有读取到数据或文件为空" << std::endl;
	}
	return 0;
}


//用来读取sp3的函数
// 输入文件名 inputFilename
int myReadSp3_C(std::string inputFilename) {
	std::string outputFilename = "sp3output_C.csv"; // 输出文件名
	int timestamp = -5;
	std::vector<SatelliteData>data;
	std::ifstream infile(inputFilename);
	if (!infile.is_open())
	{
		std::cerr << "无法打开文件：" << inputFilename << std::endl;
		data = {};
	}
	std::string line;
	bool firstStarEnded = false;
	while (std::getline(infile, line))
	{
		if (!firstStarEnded)
		{
			if (line.find("/* PCV:IGS20      OL/AL:FES2014b NONE     YN ORB:CoN CLK:CoN  ") != std::string::npos)
			{
				firstStarEnded = true;
			}
			continue;
		}
		if (line.rfind("*", 0) == 0)
		{
			timestamp += 5;
		}
		else if (line.rfind("P", 0) == 0)
		{
			//检查字符串是否包含子字符串"PC",0表示从零开始
			size_t found = line.find("PC");
			if (found != std::string::npos)
			{
				//数据行
				SatelliteData satellite;
				satellite.timeStamp = timestamp;
				std::istringstream ss(line);
				ss >> satellite.satelliteId >> satellite.x >> satellite.y >> satellite.z >> satellite.clockBias;
				data.push_back(satellite);
			}
			else
				continue;
		}
	}
	infile.close();

	if (!data.empty()) {
		saveDataToCSV(outputFilename, timestamp, data);
		std::cout << "数据已保存到 " << outputFilename << std::endl;
	}
	else {
		std::cerr << "没有读取到数据或文件为空" << std::endl;
	}
	return 0;
}

void menuDisplay()
{
	const int menuWidth = 100;
	const std::string title = " Main Menu ";
	const std::vector<std::string> options = {
		"1. Run the GPS satellite position calculation program",
		"2. Run the Beidou satellite position calculation program",
		"3. Run the SPP(single point positioning) using GPS",
		"4. Run the SPP using Iono-Free Combination",
		"5. Exit"
	};

	printBorderLine(menuWidth);
	printMenuItem(title, menuWidth);
	printMiddleBorderLine(menuWidth);

	for (const std::string& option : options)
	{
		printMenuItem(option, menuWidth);
	}

	printBottomBorderLine(menuWidth);
}

int getUserSelection()
{
	int selection;
	while (true)
	{
		std::cout << "Please enter your choice (1-5): ";
		std::cin >> selection;

		if (std::cin.fail())
		{
			std::cin.clear(); // clear the error flag
			std::cin.ignore(1024, '\n'); // discard invalid input
			std::cout << "Invalid input. Please enter a number between 1 and 4." << std::endl;
		}
		else if (selection < 1 || selection > 4)
		{
			std::cout << "Invalid choice. Please enter a number between 1 and 4." << std::endl;
		}
		else
		{
			return selection;
		}
	}
}

void handleUserSelection(int selection)
{
	switch (selection)
	{
		case 1:
		{
			std::cout << "Running the GPS satellite position calculation program..." << std::endl;
			myReadSp3("COD0MGXFIN_20240500000_01D_05M_ORB.SP3");
			vector<double> tVector = tVectorGenerate();
			vector<nav_body>g_vector = NDataRead_G("BRDC00IGS_R_20240500000_01D_MN.rnx");
			vector<NavData>nav_vectors = myReadRNX(g_vector, tVector);
			//GPS卫星总数据
			std::cout << "计算得到的卫星的数据是" << endl;
			for (NavData nav : nav_vectors)
			{
				cout << "在2月19日" << nav.time << "秒（时刻），卫星G" << nav.sPRN << "的位置是" << "x：" << nav.X << "y：" << nav.Y << "z：" << nav.Z << endl;
				cout << "此时卫星的钟差是" << nav.clockBias << endl;
				cout << "--------------------------------分割线--------------------------------" << endl;
			}
		}
			break;
		case 2:
		{
			std::cout << "Running the Beidou satellite position calculation program..." << std::endl;
			myReadSp3_C("COD0MGXFIN_20240500000_01D_05M_ORB.SP3");
			vector<double> tVector = tVectorGenerate();
			vector<nav_body>c_vector = NDataRead_C("BRDC00IGS_R_20240500000_01D_MN.rnx");
			vector<NavData>nav_vectors = myReadRNX_C(c_vector, tVector);
			//GPS卫星总数据
			std::cout << "计算得到的BDS卫星的数据是" << endl;
			for (NavData nav : nav_vectors)
			{
				cout << "在2月19日" << nav.time << "秒（时刻），卫星C" << nav.sPRN << "的位置是" << "x：" << nav.X << "y：" << nav.Y << "z：" << nav.Z << endl;
				cout << "此时卫星的钟差是" << nav.clockBias << endl;
				cout << "--------------------------------分割线--------------------------------" << endl;
			}
		}
			break;
		case 3:
		{
			std::cout << "Running the SPP(single point positioning) using GPS..." << std::endl;

			// Add the SPP calculation program logic here
			try {
							//GPS卫星总数据
			vector<nav_body> g_vector;
			//读入N文件开始计算卫星位置
			g_vector = NDataRead_G("BRDC00IGS_R_20240500000_01D_MN.rnx");
			OData data = readOFile("ACRG00GHA_R_20240500000_01D_30S_MO.rnx");

			pnav_head pnavHead = readNavDataHead("BRDC00IGS_R_20240500000_01D_MN.rnx");
			// 输出结果或进一步处理
			positioning(data, g_vector, pnavHead);

			}
			catch (const std::exception& e) {
				std::cerr << e.what() << std::endl;
			}
		}
			break;
		case 4:
		{
			std::cout << "Running the SPP(single point positioning) using Iono-Free Combination" << std::endl;
			// Add the SPP calculation program logic here
			try {
				//GPS卫星总数据
				vector<nav_body> g_vector;
				//读入N文件开始计算卫星位置
				g_vector = NDataRead_G("BRDC00IGS_R_20240500000_01D_MN.rnx");
				OData data = readOFile("ACRG00GHA_R_20240500000_01D_30S_MO.rnx");
				pnav_head pnavHead = readNavDataHead("BRDC00IGS_R_20240500000_01D_MN.rnx");
				cout << "执行成功" << endl;
				// 输出结果或进一步处理
				positioning_iono_free(data, g_vector, pnavHead);

			}
			catch (const std::exception& e) {
				std::cerr << e.what() << std::endl;
			}
		}
		break;
		case 5:
			std::cout << "Exiting the program..." << std::endl;
			exit(0);
		default:
			std::cerr << "Error: Invalid selection" << std::endl;
	}
}

void runMyCode()
{
	menuDisplay();
	int selection = getUserSelection();
	handleUserSelection(selection);
}

int main()
{
	runMyCode();
	return 0;
}


