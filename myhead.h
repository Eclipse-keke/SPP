#pragma once

#include<vector>
#include<chrono>
#include<iostream>
#include<math.h>
#include<cmath>
#include<chrono>
#include <sstream>
#include <string>
#include <iomanip>  
#include<fstream> //ifstream ����Ҫ�Ŀ�
#include<cctype>//����ָʾ::ispace����
#include <cstdio>
#include<stdio.h>
#include<cmath>
#include<algorithm>
#include <stdexcept>
#ifdef _WIN32
#include<windows.h>
#else
#include<unistd.h>
#endif

//The product of the gravitational constant G and the total mass M of the Earth
#define GM 3.986005e14
//The speed of light
#define C 299792458.0
//rotational angular velocity of the earth
#define OMEGA_E 7.2921151467e-5 //unit:rad/s
#define MIU 3.986004418e14 //��λ:m3/s2 ��BDCS����ϵ�µĵ���������������ʵ����GM
#define OMEGA_SELF_ROTATE 7.2921150e-5 //BDCS����ϵ�µĵ�����ת���ٶ�
#define PI 3.1415926535897932384626
#define MAXRINEX 84 //һ���������������80�Ϳ��ԣ�
#define  F  -4.442807633e-10

enum Mode {
	NO_CORRECTION = 1,
	IONO_CORRECTION = 2,
	TROPO_CORRECTION = 3,
	TROPO_IONO_CORRECTION = 4
};

typedef struct nav_head
{
	double ver;//rinex �汾��
	char type[16];//��ȡ����������
	double ION_alpha[4];//8����������
	double ION_beta[4];
}nav_head, * pnav_head;

//����N�ļ����ݽṹ��
typedef struct nav_body
{
	//The monument of Satillite, can be 'P':GPS, 'C':Beidou
	char monument;
	//The content of the first line of the data block
	int sPRN; //PRN of a satellite
	//��GNSS��ȫ�򵼺�����ϵͳ���У�
	// PRN����α�����������Pseudo-Random Noise��
	// ����д��PRN��һ�������α��������У�
	//���ڱ�ʶ�����ֲ�ͬ�����ǣ���ÿ�����ǵ�ID��
	int TOC_Y;//Year
	int TOC_M;//Month
	int TOC_D;//Day
	int TOC_H;//Hour
	int TOC_Min;//Minute
	int TOC_Sec;//Second
	double sa0;//The Clock Error of the Satellite, �����Ӳ�
	double sa1;//The Satellite clock deviation������ƫ
	double sa2;//The Satellite clock drift������Ư


	//The Content of the 2nd line of the data block
	double IODE;//Data and Ephemeris�������� release time
	double Crs;//����뾶�����ҵ��͸�������������λ��m��
	double delta_n;//����ƽ���˶����������ֵ�Ĳ�ֵ����λ��rad/s����n�����˹�������صĳ����㶯
	double M0;//�ο�ʱ���ƽ�����


	//The Content of the 3rd line of the data block
	double Cuc;//ά�ȷ��ǵ����ҵ��͸�����������Rad��
	double e;//���ƫ����
	double Cus;//������ǵ����ҵ��͸�����������Rad��
	double sqrtA;//������ƽ����

	//The Content of the 4th line of the Data block
	double TOE;//�����Ĳο�ʱ�̣�GPS�����룩
	double Cic;// �����ǵ����ҵ��͸�����������rad��
	double OMEGA;//�ο�ʱ�̵�������ྭ
	double Cis; //ά����ǵ����ҵ��͸�����������rad��


	//���ݿ���������ݣ�
	double i0;//�ο�ʱ��Ĺ�����(rad)
	double Crc;//���ƽ�������ҵ��͸���������(m)
	double omega;//���ص�Ǿ�
	double deltaomega;//������ྭ�仯��(rad)

	//���ݿ���������ݣ�
	double IDOT;//i�ı仯��
	double L2code;//L2�ϵ���
	double GPSweek;//GPS��,��TOEһͬ��ʾ
	double L2Pflag;//L2,p�����ݱ��

	//���ݿ����������
	double sACC;//���Ǿ���
	double sHEA;//���ǽ���״̬
	double TGD;//Time Group Delay�ź��������ڵ�Ⱥ�Ӳ�
	double IODC;//�ӵ���������

	//���ݿ�ڰ�������
	double TTN;//���ķ���ʱ��
	double fit;//�������
	double spare1;//��
	double spare2;//��
}nav_body, * pnav_body;




//����Ҫ����csv������
// ���� nav_b �ṹ����Ӧ�������Ѿ�����
struct NavData {
	std::string monument;
	int sPRN=0;
	double time = 0;
	double tk =0;
	double X =0;
	double Y =0;
	double Z =0;
	double clockBias = 0;
};

//����sp3�ļ��Ľṹ�嶨��
struct SatelliteData {
	int timeStamp = 0;
	std::string satelliteId;
	double x = 0;
	double y = 0;
	double z = 0;
	double clockBias = 0;
};


//���ջ�λ�õĽṹ��
struct rcvposCoor
{
	double x;
	double y;
	double z;
};

//O�ļ����ڵĽṹ�壬�����ꡢ�¡��ա�ʱ���֡���
struct ODate
{
	int year;
	int month;
	int day;
	int hour;
	int minute;
	double second;
};


//����O�ļ������ݴ洢
struct OData
{
	ODate date;//���ڣ�
	std::vector<double>epoch;//��Ԫ
	std::vector<std::string>type;//
	std::string station;//վ��
	std::vector<std::vector<double>>data;
	std::vector<int> index;
	std::vector<rcvposCoor>rcvpos;//receive position�����ջ�λ��
	std::vector<std::vector<double>>time;
};


//SPP��ͷ�ļ�
struct ECEF
{
	double x;
	double y;
	double z;
};

struct BLH
{
	double B;
	double L;
	double H;
};

struct ENU
{
	double E;
	double N;
	double U;
};

//rΪ�����򾶣�AΪ���Ƿ�λ�ǣ�hΪ���ǵĸ߶Ƚ�
struct RAH
{
	double R;
	double A;
	double H;
};

// ���ڴ洢��ά�����Ľṹ
struct Vector3D {
	double x;
	double y;
	double z;
};

// ���ڴ洢��ά�����Ľṹ
struct Vector2D {
	double x;
	double y;
};

//���ڴ洢4ά�����ṹ
struct Vector4D
{
	double x;
	double y;
	double z;
	double br;
};
// ���ڴ洢λ�ý��
struct PositioningResult 
{
	std::vector<Vector3D> blh; // ����ʹ�� {latitude, longitude, height} ��ʾ
	std::vector<Vector3D> xyz; // ECEF �µ�����
	std::vector<std::vector<double>> bs; // �����ӲSatellite Clock Bias
	std::vector<double> br; // ���ջ��Ӳ�
	std::vector<std::vector<double>> elevation; // elevationAngle
};

// ���ڴ洢�������
struct DistError {
	std::vector<double> horizontal; // ��λ��Horizontal
	std::vector<double> EW; // ��λ����������
	std::vector<double> NS; // ��λ���ϱ�����
	std::vector<double> height; // ��λ���̷߳���
};

// ���ڴ洢����ģ��
struct Model {
	std::vector<std::vector<double>> tropo; // �������ӳ�
	std::vector<std::vector<double>> iono; // ������ӳ�
	std::vector<double> hdop; // �߳̾���˥������ HDOP
};

// ���ڴ洢 N �ļ���������Ҳ�������ǵ�λ�ã�
struct NavDataResult {
	std::vector<Vector3D> eph; // ����
	std::vector<int> index; // ���
	std::vector<double> inoprm; // ��������
};

/*main.cpp*/
int getrow(FILE* fp_nav);
std::vector<nav_body> NDataRead_G(const char* fileName);
NavData myComputeLocationOfSatellite(nav_body nav_b, double tk);
void readrinex_head(FILE* fp_nav, pnav_head nav_h);
int readrinex_body(const char* fileName, std::vector<nav_body>& nav_b_vector,int n_n,char choose);
void dataBlocksRead(std::vector<std::vector<std::string>>dataBlocks, std::vector<nav_body>& nav_b_vector,char choose);
nav_body findTrueOne(double t, std::vector<nav_body> vessel);
std::vector<NavData> myReadRNX(std::vector<nav_body>g_vector, std::vector<double>tVector);
std::vector<NavData> myReadRNX_C(std::vector<nav_body>c_vector, std::vector<double>tVector);
static double strtonum(const std::string& buff, int i, int n);
static double strtonum(const char* buff, int i, int n);

/*compute tools*/
double toc_compute(double toe, double TOC_H, double TOC_Min, double TOC_Sec);
double computeAverageAngularVelocity(double squareRootA, double delta_n);
double computeInstantMeanAnomaly(double M0, double n, double tk);
double computeEccentricAnomaly(double M, double e);
double solveTrueAnomaly(double E, double e);
double solveArgumentOfLatitude(double omega, double f);
void solvePerturbationCorrection(double Cuc, double Cus, double Crc, double Crs, double Cic, double Cis, double* delta_u, double* delta_r, double* delta_i, double u_pie);
void perturbationCorrection(double sqrtA, double* u, double* r, double* i, double IDOT, double E, double tk, double e, double u_pie, double i0, double delta_u, double delta_r, double delta_i);
void computexyInOrbitPlane(double r, double u, double* x, double* y);
double computeMomentOMEGA(double tk, double deltaomega, double OMEGA, double TOE);
void solveSatillitePosition(double L, double i, double x, double y, double* X, double* Y, double* Z);
double clockBiasCompute(nav_body nav_b, double tk);
double norm(const std::vector<double>& vec);
RAH sateliite_rah(ECEF receiver, ECEF satellite);
ENU ecef2enu(ECEF receiver, ECEF satellite);
BLH xyz2blh(ECEF ecef);
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix);
std::vector<std::vector<double>> createIdentityMatrix(int n);
std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);

/*sp3Read.cpp*/
std::vector<SatelliteData> readDataFromFile(const std::string& filename, int& timestamp);
void saveDataToCSV(const std::string& filename, int&timestamp, const std::vector<SatelliteData>& data);

/*readOFiles.cpp*/
OData readOFile(const std::string& r_o_name);

/*SPP.cpp*/
void positioning(const OData& obs, std::vector<nav_body> g_vector, pnav_head pnavHead);
void removeOutliers(std::vector<double>& cP1, std::vector<NavData>& satpos_m, const rcvposCoor& refpos);
std::vector<double> computeRefDistance(const rcvposCoor& refpos, const std::vector<NavData>& satpos_m);
void positioning_iono_free(const OData& obs, std::vector<nav_body> g_vector, pnav_head pnavHead);

/*tropoAndIono.cpp*/
double klobuchar_model(double fi, double lambda, double elev, double tow, const double ionprm[2][4]);
void EstimateTropDalay(double Latitude, double Height, double DOY, double& d_hyd, double& d_wet);


/*beautifulMenu.cpp*/
void printCentered(const std::string& text, int width);
void printBorderLine(int length);
void printMiddleBorderLine(int length);
void printBottomBorderLine(int length);
void printMenuItem(const std::string& item, int width);