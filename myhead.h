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
#include<fstream> //ifstream 类需要的库
#include<cctype>//用来指示::ispace函数
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
#define MIU 3.986004418e14 //单位:m3/s2 是BDCS坐标系下的地心引力常数，其实就是GM
#define OMEGA_SELF_ROTATE 7.2921150e-5 //BDCS坐标系下的地球自转角速度
#define PI 3.1415926535897932384626
#define MAXRINEX 84 //一行最大字数（大于80就可以）
#define  F  -4.442807633e-10

enum Mode {
	NO_CORRECTION = 1,
	IONO_CORRECTION = 2,
	TROPO_CORRECTION = 3,
	TROPO_IONO_CORRECTION = 4
};

typedef struct nav_head
{
	double ver;//rinex 版本号
	char type[16];//读取的数据类型
	double ION_alpha[4];//8个电离层参数
	double ION_beta[4];
}nav_head, * pnav_head;

//创建N文件数据结构体
typedef struct nav_body
{
	//The monument of Satillite, can be 'P':GPS, 'C':Beidou
	char monument;
	//The content of the first line of the data block
	int sPRN; //PRN of a satellite
	//在GNSS（全球导航卫星系统）中，
	// PRN代表“伪随机噪声”（Pseudo-Random Noise）
	// 的缩写。PRN是一种特殊的伪随机码序列，
	//用于标识和区分不同的卫星，即每颗卫星的ID。
	int TOC_Y;//Year
	int TOC_M;//Month
	int TOC_D;//Day
	int TOC_H;//Hour
	int TOC_Min;//Minute
	int TOC_Sec;//Second
	double sa0;//The Clock Error of the Satellite, 卫星钟差
	double sa1;//The Satellite clock deviation卫星钟偏
	double sa2;//The Satellite clock drift卫星钟漂


	//The Content of the 2nd line of the data block
	double IODE;//Data and Ephemeris（星历） release time
	double Crs;//轨道半径的正弦调和改正项的振幅（单位：m）
	double delta_n;//卫星平均运动速率与计算值的差值（单位：rad/s）△n包括了轨道参数ω的长期摄动
	double M0;//参考时间的平近点角


	//The Content of the 3rd line of the data block
	double Cuc;//维度幅角的余弦调和改正项的振幅（Rad）
	double e;//轨道偏心率
	double Cus;//轨道幅角的正弦调和改正项的振幅（Rad）
	double sqrtA;//长半轴平方根

	//The Content of the 4th line of the Data block
	double TOE;//星历的参考时刻（GPS周内秒）
	double Cic;// 轨道倾角的余弦调和改正项的振幅（rad）
	double OMEGA;//参考时刻的升交点赤经
	double Cis; //维度倾角的正弦调和改正项的振幅（rad）


	//数据块第五行内容：
	double i0;//参考时间的轨道倾角(rad)
	double Crc;//轨道平径的余弦调和改正项的振幅(m)
	double omega;//近地点角距
	double deltaomega;//升交点赤经变化率(rad)

	//数据块第六行内容：
	double IDOT;//i的变化率
	double L2code;//L2上的码
	double GPSweek;//GPS周,与TOE一同表示
	double L2Pflag;//L2,p码数据标记

	//数据块第七行内容
	double sACC;//卫星精度
	double sHEA;//卫星健康状态
	double TGD;//Time Group Delay信号在卫星内的群延差
	double IODC;//钟的数据龄期

	//数据块第八行内容
	double TTN;//电文发送时间
	double fit;//拟合区间
	double spare1;//空
	double spare2;//空
}nav_body, * pnav_body;




//这是要导出csv表格里的
// 假设 nav_b 结构和相应的数据已经定义
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

//这是sp3文件的结构体定义
struct SatelliteData {
	int timeStamp = 0;
	std::string satelliteId;
	double x = 0;
	double y = 0;
	double z = 0;
	double clockBias = 0;
};


//接收机位置的结构体
struct rcvposCoor
{
	double x;
	double y;
	double z;
};

//O文件日期的结构体，包括年、月、日、时、分、秒
struct ODate
{
	int year;
	int month;
	int day;
	int hour;
	int minute;
	double second;
};


//这是O文件的数据存储
struct OData
{
	ODate date;//日期：
	std::vector<double>epoch;//历元
	std::vector<std::string>type;//
	std::string station;//站名
	std::vector<std::vector<double>>data;
	std::vector<int> index;
	std::vector<rcvposCoor>rcvpos;//receive position：接收机位置
	std::vector<std::vector<double>>time;
};


//SPP的头文件
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

//r为卫星向径，A为卫星方位角，h为卫星的高度角
struct RAH
{
	double R;
	double A;
	double H;
};

// 用于存储三维向量的结构
struct Vector3D {
	double x;
	double y;
	double z;
};

// 用于存储二维向量的结构
struct Vector2D {
	double x;
	double y;
};

//用于存储4维向量结构
struct Vector4D
{
	double x;
	double y;
	double z;
	double br;
};
// 用于存储位置结果
struct PositioningResult 
{
	std::vector<Vector3D> blh; // 可以使用 {latitude, longitude, height} 表示
	std::vector<Vector3D> xyz; // ECEF 下的坐标
	std::vector<std::vector<double>> bs; // 卫星钟差：Satellite Clock Bias
	std::vector<double> br; // 接收机钟差
	std::vector<std::vector<double>> elevation; // elevationAngle
};

// 用于存储距离误差
struct DistError {
	std::vector<double> horizontal; // 定位误差：Horizontal
	std::vector<double> EW; // 定位误差：东西方向
	std::vector<double> NS; // 定位误差：南北方向
	std::vector<double> height; // 定位误差：高程方向
};

// 用于存储改正模型
struct Model {
	std::vector<std::vector<double>> tropo; // 对流层延迟
	std::vector<std::vector<double>> iono; // 电离层延迟
	std::vector<double> hdop; // 高程精度衰减因子 HDOP
};

// 用于存储 N 文件计算结果（也就是卫星的位置）
struct NavDataResult {
	std::vector<Vector3D> eph; // 星历
	std::vector<int> index; // 编号
	std::vector<double> inoprm; // 电离层参数
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