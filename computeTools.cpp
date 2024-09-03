#define _CRT_SECURE_NO_WARNINGS
#include"myhead.h"

/*
This tool uses broadcast ephemeris to compute the specific location of satellite
*/



/*
计算时间t
d240101//2024年1月1日，是星期一
tk 从参考历元开始算的时间, tk = t - toe, t为观测时刻
tyear
tmonth
tday
thour
tsec
*/
double toc_compute(double toe,double TOC_H, double TOC_Min, double TOC_Sec)
{
    double d240219 = 1;//因为是星期一所以是1
    //double tyear = TOC_Y * 86400 * 30.6001 * 362.25;
   // double tmonth = TOC_M * 86400 * 30.6001;
  //  double tday = TOC_D * 86400;
    double thour = TOC_H * 3600;
    double tmin = TOC_Min * 60;
    double tsec = TOC_Sec;
    //toc就是用来计算钟差的参考时刻
    double toc = d240219 * 86400 + thour + tmin + tsec;
    return toc;
}


/*
This is for computing the average angular velocity of satellite.
Input:
第三行
double sqrtA;//长半轴平方根
第四行
double TOE;//星历的参考时刻（GPS周内秒）
第2行
double delta_n;//卫星平均运动速率与计算值的差值（单位：rad/s）△n包括了轨道参数ω的长期摄动
*/
double computeAverageAngularVelocity(double squareRootA, double delta_n)
{
	//卫星运动的平均角速度
	double n = 0;
	//参考时刻Toe的平均角速度
	double n0 = 0;
	n0 = std::sqrt(GM) / (std::pow(squareRootA, 3));
	n = n0 + delta_n;
	return n;
}


/*
* This function is for computing the mean anomaly of satellite at the moment of observation.
* 	double M0;//参考时间的平近点角
* input: M0, Toe, n,t
* t is the moment of the observation
*/
double computeInstantMeanAnomaly(double M0, double n,double tk)
{
	double M = 0;
	M = M0 + n * tk;
	return M;
}




/*
在椭圆的参数方程x=acosθ ， y=bsinθ中，参数角θ即为偏近点角。
Using rad to compute
input:
e and M
使用迭代解法，E的初值可以是M，迭代精度1e-12，最大迭代次数1000
*/
double computeEccentricAnomaly(double M,double e)
{
    double E_old = M;
    double E = M;
    double dE = 1.0;
    int maxIteration = 1000;
    int iteration = 0;
    while (dE > 1e-12&&iteration<=1000)
    {
        E = M + e * sin(E_old);
        dE = abs(E - E_old);
        E_old = E;
        iteration++;
    }
    return E;
}



/*
compute True Anomaly
input: e
*/
double solveTrueAnomaly(double E, double e)
{
    double f = 0;
    f = 2.0 * atan(    sqrt(1 - e * e) * tan(E / 2.0)      /      (1 - e)    );
    //f = std::atan(std::sqrt(1-e*e)*std::sin(E)/(std::cos(E) - e));
    return f;
}

/*
solve the argument of latitude
input: omega and f
double omega;//近地点角距,第五行
*/
double solveArgumentOfLatitude(double omega, double f)
{
    double u_pie = omega + f;
    return u_pie;
}


/*
solve the perturbation correction delta_u, delta_r,delta_i
input:
Cuc, Cus, Crc, Crs, Cic,Cis
*/
void solvePerturbationCorrection(double Cuc, double Cus, double Crc, double Crs, double Cic, double Cis, double* delta_u, double*delta_r,double*delta_i,double u_pie)
{
    *delta_u = Cuc * std::cos(2 * u_pie) + Cus * std::sin(2 * u_pie);
    *delta_r = Crc * std::cos(2 * u_pie) + Crs * std::sin(2 * u_pie);
    *delta_i = Cic * std::cos(2 * u_pie) + Cis * std::sin(2 * u_pie);
    return;
}


/*
Perturbation correction for u_pie, r_pie. i0
input:
 数据块第六行内容：
 double IDOT;//i的变化率
*/
void perturbationCorrection(double sqrtA, double* u, double* r, double* i,double IDOT, double E,double tk, double e,double u_pie, double i0,double delta_u, double delta_r, double delta_i)
{
    double a = sqrtA * sqrtA;
    *u = u_pie + delta_u;
    *r = a * (1 - e * std::cos(E)) + delta_r;
    *i = i0 + delta_i + IDOT * tk;
    return;
}




/*
solve the position of satillite in coordinates of satillite orbit plane.
*/
void computexyInOrbitPlane(double r, double u,double *x, double *y)
{
    *x = r * std::cos(u);
    *y = r * std::sin(u);
}

/*
solve the longitude L of ascending node at the observation moment
using:
    	double TOE;//星历的参考时刻（GPS周内秒）
    double deltaomega;//升交点赤经变化率(rad)
        double OMEGA;
        //注意：广播星历中给出的并不是参考时刻的toe时的升交点赤经，而是该值与本周起始时刻的格林尼治恒星时GAST_week之差
        也就是
        OMEGA = OMEGA_toe- GAST_week

        t在这里为本周时间
*/
double computeMomentOMEGA(double tk, double deltaomega, double OMEGA, double TOE)
{
    //计算观测瞬间的升交点赤经
    double L = OMEGA + (deltaomega - OMEGA_E) * tk - OMEGA_E * TOE;
    return L;
}


/*
计算卫星在瞬时地球坐标系中的位置
输入：升交点的大地经度L，轨道平面的倾角i，卫星在轨道平面的位置x,y
输出：X Y Z
*/
void solveSatillitePosition(double L, double i,double x, double y,double *X, double *Y, double *Z)
{
    *X = (x * std::cos(L) - y * std::cos(i) * std::sin(L));
    *Y =(x * std::sin(L) + y * std::cos(i) * std::cos(L));
    *Z = (y * std::sin(i));
}



/*
由广播星历计算卫星钟差
*/
double clockBiasCompute(nav_body nav_b, double tk)
{
    //钟差计算

    //1.相对改正（tr是GPS卫星非圆形轨道而引起的相对论效应的修正项）
    //1.1计算Eccentric Anomaly
    double n = computeAverageAngularVelocity(nav_b.sqrtA, nav_b.delta_n);
    double M = computeInstantMeanAnomaly(nav_b.M0, n, tk);
    double EA = computeEccentricAnomaly(M, nav_b.e);

    double r_c = F * nav_b.e * nav_b.sqrtA * sin(EA);

    //1.2计算SV clock correction
    double t_sv = nav_b.sa0 + nav_b.sa1 * tk + nav_b.sa2 * tk * tk + r_c;

    //1,3得到卫星钟差
    double satclock = t_sv - nav_b.TGD;

    return satclock;
}


/*
计算向量的二范数
*/
double norm(const std::vector<double>& vec)
{
    double sum = 0.0;
    for (const auto& value : vec)
    {
        sum += value * value;
    }
    return std::sqrt(sum);
}



//ECEF(x,y,z) To BLH(B,L,H)
BLH xyz2blh(ECEF ecef)
{
    //以下是椭球参数
    double a = 6378137.0;
    double f = 1.0 / 298.257223563;
    double e2 = f * (2 - f);
    BLH blh = {};
    double R0 = sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    double R1 = sqrt(ecef.x * ecef.x + ecef.y * ecef.y + ecef.z * ecef.z);
    //solve the Longitude
    blh.L = atan2(ecef.y, ecef.x);
    //Using iterations to solve the latitude of geodesy and Height of geodesy
    double N = a;
    double H = R1 - a;
    double B = atan2(ecef.z * (N + H), R0 * (N * (1 - e2) + H));

    //进行迭代，这是每次迭代更新的量
    double dH;
    double dB;
    do
    {
        dH = N;
        dB = B;
        N = a / sqrt(1 - e2 * sin(B) * sin(B));
        H = R0 / cos(B) - N;
        B = atan2(ecef.z * (N + H), R0 * (N * (1 - e2) + H));
    } while (fabs(H - dH) > 1e-3 && fabs(B - dB) > 1e-9);

    blh.B = B;
    blh.H = H;
    return blh;
}


//ecef(xyz) To 东北天坐标系ENU
//输入：测站在ecef的坐标xr,yr,zr, 卫星在ecef的坐标xs,ys,zs
ENU ecef2enu(ECEF receiver, ECEF satellite)
{
    ENU enu = {};
    BLH blh = xyz2blh(receiver);
    double sinL = sin(blh.L);
    double cosL = cos(blh.L);
    double sinB = sin(blh.B);
    double cosB = cos(blh.B);
    double dx = satellite.x - receiver.x;
    double dy = satellite.y - receiver.y;
    double dz = satellite.z - receiver.z;

    enu.E = -sinL * dx + cosL * dy;
    enu.N = -sinB * cosL * dx - sinB * sinL * dy + cosB * dz;
    enu.U = cosB * cosL * dx + cosB * sinL * dy + sinB * dz;
    return enu;
}

//求卫星高度角
//A：方位角
//r为卫星向径，A为卫星方位角，h为卫星的高度角
RAH sateliite_rah(ECEF receiver, ECEF satellite)
{
    RAH rah;
    BLH user_blh = xyz2blh(receiver);
    BLH sat_blh = xyz2blh(satellite);
    double Lat = user_blh.B;
    double Lon = user_blh.L;
    double Lat_s = sat_blh.B;
    double Lon_s = sat_blh.L;
    double R[3][3] = {
    {-sin(Lon),                cos(Lon),                0},
    {-sin(Lat) * cos(Lon), -sin(Lat) * sin(Lon),  cos(Lat)},
    {cos(Lat) * cos(Lon),   cos(Lat) * sin(Lon),  sin(Lat)}
    };

    double Rs[3] = {
                satellite.x - receiver.x,
        satellite.y - receiver.y,
        satellite.z - receiver.z
    };
    double RL[3];
    for (int i = 0; i < 3; ++i) {
        RL[i] = 0;
        for (int j = 0; j < 3; ++j) {
            RL[i] += R[i][j] * Rs[j];
        }
    }

    double Xl = RL[0];
    double Yl = RL[1];
    double Zl = RL[2];

    rah.H = atan2(Zl, sqrt(Xl * Xl + Yl * Yl)) * 180 / PI;
    rah.A = atan2(cos(Lat_s) * sin(Lon_s - Lon),
        cos(Lat) * sin(Lat_s) - sin(Lat) * cos(Lat_s) * cos(Lon_s - Lon));
    rah.A = rah.A * 180 / PI;
    rah.R = 0;
    return rah;
}


/*矩阵运算*/
// 矩阵的转置
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) 
{
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows, 0.0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

//矩阵的乘法
std::vector<std::vector<double>> matrixMultiply(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B)
{
    int rowsA = A.size();
    int colsA = A[0].size();
    int rowsB = B.size();
    int colsB = B[0].size();
    std::vector<std::vector<double>> result(rowsA, std::vector<double>(colsB, 0.0));

    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// 辅助函数：创建单位矩阵
std::vector<std::vector<double>> createIdentityMatrix(int n)
{
    std::vector<std::vector<double>> identity(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }
    return identity;
}

// 高斯-约当消元法求逆
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix)
{
    int n = matrix.size();
    std::vector<std::vector<double>> result = createIdentityMatrix(n);

    // 复制输入矩阵，以免修改原始矩阵
    std::vector<std::vector<double>> A = matrix;

    // 高斯-约当消元法
    for (int i = 0; i < n; ++i) {
        // 将对角线元素归一
        double divisor = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= divisor;
            result[i][j] /= divisor;
        }

        // 将除了对角线元素外的其它列元素置零
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                    result[k][j] -= factor * result[i][j];
                }
            }
        }
    }
    return result;
}




/*
计算卫星在协议地球坐标系中的位置
*/

