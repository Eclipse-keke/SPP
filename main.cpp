#define _CRT_SECURE_NO_WARNINGS
#include"myhead.h"
using namespace std;


void saveToCSV(const std::vector<NavData>& navData,const std::string& filename)
{
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "�޷����ļ�: " << filename << std::endl;
		return;
	}

	// д��CSV��ͷ
	file << "���Ǳ�ʶ,PRN,tk,X,Y,Z,clockBias" << std::endl;

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
����t�Ļ���
*/
vector<double> tVectorGenerate()
{
	vector<double> tVector;
	double d240219 = 1;//��Ϊ������һ������1
	for (int i = 0; i <= 1440; i += 5)
	{
		double t = d240219 * 86400 + i*60;
		tVector.push_back(t);
	}
	return tVector;
}


/*
���ܣ��ҵ�t-toc��С����һ������
����t����tVector�����ģ����������Ǳ�ţ������Ǹ��������������ݵ�����
����Ǹ���Ӧ�����������Ǳ�ţ�ʱ�����
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
		nav_b_indexes.push_back(i); // ��¼ nav_b ������
		if (delta_t < delta_t_min) // �ҵ���С��ֵ
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
	nav_h = (pnav_head)malloc(sizeof(nav_head));//��N�ļ�ͷ���ٿռ�
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




//�ļ���ȡ������
//�������ļ����� / ���·��
//��ȡN�ļ���GPS�������ݣ����õĽ���Ǻܶ��nav bodyȻ��洢��һ��vector��
vector<nav_body> NDataRead_G(const char* fileName)
{
	//Read Head File
	FILE* fp_nav = NULL;//���������ļ�ָ��
	pnav_body nav_b = NULL;

	//N�ļ���ʼ��ȡ
	fp_nav = fopen(fileName, "r");
	int n = getrow(fp_nav); //��ȡ���������ļ�������
	std::cout << "���������ļ���������" << endl;
	cout << n << endl;
	rewind(fp_nav);//Move the file pointer to the starting position

	nav_b = (pnav_body)malloc(sizeof(nav_body));//malloc space for body of File N.

	fclose(fp_nav);

	//GPS����������
	vector<nav_body> g_vector;
	//Read Body File.
	readrinex_body(fileName, g_vector, n,'G');
	return g_vector;
}

//�������ļ����� / ���·��
//��ȡN�ļ���Beidou�������ݣ����õĽ���Ǻܶ��nav bodyȻ��洢��һ��c_vector��
vector<nav_body> NDataRead_C(const char* fileName)
{
	//Read Head File
	FILE* fp_nav = NULL;//���������ļ�ָ��
	pnav_body nav_b = NULL;

	//N�ļ���ʼ��ȡ
	fp_nav = fopen(fileName, "r");
	int n = getrow(fp_nav); //��ȡ���������ļ�������
	std::cout << "���������ļ���������" << endl;
	cout << n << endl;
	rewind(fp_nav);//Move the file pointer to the starting position

	nav_b = (pnav_body)malloc(sizeof(nav_body));//malloc space for body of File N.

	fclose(fp_nav);

	//Beidou����������
	vector<nav_body> c_vector;
	readrinex_body(fileName, c_vector, n,'C');
	return c_vector;
}



/*
��������
*/
vector<NavData> myReadRNX(vector<nav_body>g_vector,vector<double>tVector)
{
	//GPS���ǵ����������ݣ���32���洢
	vector<vector<nav_body>>  gpsSats;
	gpsSats.resize(32);


	cout << "���" << endl;
	//��ʱ��ת��Ϊgpsʱ
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
			//1. ��tVector���ÿһ��t��Ѱ�������t�����toe��Ӧ��һ�������Ҳ���ǣ�t-toe>����t-toe��С
			trueParams = findTrueOne(t, gpsSats[i]);
			//2.�����t���㣬����t-toe
			NavData navData;
			tk = t - trueParams.TOE;
			navData = myComputeLocationOfSatellite(trueParams, tk);//t�������룬һ�����ڵĴ�����һ86400*7
			navData.time = t-86400;
			data.push_back(navData);
		}
	}
	saveToCSV(data, "brdcResult.csv");
	std::cout << "�����ѱ��浽 brdcResult.csv" << std::endl;
	return data;
}

/*
��������2��For Beidou
*/
vector<NavData> myReadRNX_C(vector<nav_body>c_vector, vector<double>tVector)
{
	//�������ǵ����������ݣ���62���洢
	vector<vector<nav_body>>  bdsSats;
	bdsSats.resize(62);


	cout << "���" << endl;
	//��ʱ��ת��Ϊgpsʱ
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
			//1. ��tVector���ÿһ��t��Ѱ�������t�����toe��Ӧ��һ�������Ҳ���ǣ�t-toe>����t-toe��С
			trueParams = findTrueOne(t, bdsSats[i]);
			//2.�����t���㣬����t-toe
			NavData navData;
			tk = t - trueParams.TOE-14;//��������Ҫ��ȥ14
			navData = myComputeLocationOfSatellite(trueParams, tk);//t�������룬һ�����ڵĴ�����һ86400*7
			navData.time = t - 86400;
			data.push_back(navData);
		}
	}
	saveToCSV(data, "brdcResult_C.csv");
	std::cout << "�����ѱ��浽 brdcResult_C.csv" << std::endl;
	return data;
}


/*
���Ǿ����������㺯��
���룺���ǲ����ṹ�壬ʱ��t
�����������˲ʱ�ع�����ϵ�µ�����
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
	cout << "������ʱ����"<<nav_b.TOC_Y<<" "<<nav_b.TOC_M<<" "<<nav_b.TOC_D<<" "<<nav_b.TOC_H<<" "<<nav_b.TOC_Min<<" "<<nav_b.TOC_Sec << endl;
	cout << "���Ǳ�ʶ:" << nav_b.monument << endl;
	cout << "PRN��" << nav_b.sPRN << endl;
	cout << "����õ�������������" << endl;
	cout << "X:"<<std::setprecision(9) << X << endl;
	cout << "Y:" << std::setprecision(9) << Y << endl;
	cout << "Z:" << std::setprecision(9) << Z << endl;
	cout << "�Ӳ���" << endl;
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

//������ȡsp3�ĺ���
// �����ļ��� inputFilename
int myReadSp3(std::string inputFilename) {
	std::string outputFilename = "sp3output.csv"; // ����ļ���
	int timestamp = -5;
	std::vector<SatelliteData> data = readDataFromFile(inputFilename, timestamp);
	if (!data.empty()) {
		saveDataToCSV(outputFilename, timestamp, data);
		std::cout << "�����ѱ��浽 " << outputFilename << std::endl;
	}
	else {
		std::cerr << "û�ж�ȡ�����ݻ��ļ�Ϊ��" << std::endl;
	}
	return 0;
}


//������ȡsp3�ĺ���
// �����ļ��� inputFilename
int myReadSp3_C(std::string inputFilename) {
	std::string outputFilename = "sp3output_C.csv"; // ����ļ���
	int timestamp = -5;
	std::vector<SatelliteData>data;
	std::ifstream infile(inputFilename);
	if (!infile.is_open())
	{
		std::cerr << "�޷����ļ���" << inputFilename << std::endl;
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
			//����ַ����Ƿ�������ַ���"PC",0��ʾ���㿪ʼ
			size_t found = line.find("PC");
			if (found != std::string::npos)
			{
				//������
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
		std::cout << "�����ѱ��浽 " << outputFilename << std::endl;
	}
	else {
		std::cerr << "û�ж�ȡ�����ݻ��ļ�Ϊ��" << std::endl;
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
			//GPS����������
			std::cout << "����õ������ǵ�������" << endl;
			for (NavData nav : nav_vectors)
			{
				cout << "��2��19��" << nav.time << "�루ʱ�̣�������G" << nav.sPRN << "��λ����" << "x��" << nav.X << "y��" << nav.Y << "z��" << nav.Z << endl;
				cout << "��ʱ���ǵ��Ӳ���" << nav.clockBias << endl;
				cout << "--------------------------------�ָ���--------------------------------" << endl;
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
			//GPS����������
			std::cout << "����õ���BDS���ǵ�������" << endl;
			for (NavData nav : nav_vectors)
			{
				cout << "��2��19��" << nav.time << "�루ʱ�̣�������C" << nav.sPRN << "��λ����" << "x��" << nav.X << "y��" << nav.Y << "z��" << nav.Z << endl;
				cout << "��ʱ���ǵ��Ӳ���" << nav.clockBias << endl;
				cout << "--------------------------------�ָ���--------------------------------" << endl;
			}
		}
			break;
		case 3:
		{
			std::cout << "Running the SPP(single point positioning) using GPS..." << std::endl;

			// Add the SPP calculation program logic here
			try {
							//GPS����������
			vector<nav_body> g_vector;
			//����N�ļ���ʼ��������λ��
			g_vector = NDataRead_G("BRDC00IGS_R_20240500000_01D_MN.rnx");
			OData data = readOFile("ACRG00GHA_R_20240500000_01D_30S_MO.rnx");

			pnav_head pnavHead = readNavDataHead("BRDC00IGS_R_20240500000_01D_MN.rnx");
			// ���������һ������
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
				//GPS����������
				vector<nav_body> g_vector;
				//����N�ļ���ʼ��������λ��
				g_vector = NDataRead_G("BRDC00IGS_R_20240500000_01D_MN.rnx");
				OData data = readOFile("ACRG00GHA_R_20240500000_01D_30S_MO.rnx");
				pnav_head pnavHead = readNavDataHead("BRDC00IGS_R_20240500000_01D_MN.rnx");
				cout << "ִ�гɹ�" << endl;
				// ���������һ������
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


