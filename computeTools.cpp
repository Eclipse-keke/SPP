#define _CRT_SECURE_NO_WARNINGS
#include"myhead.h"

/*
This tool uses broadcast ephemeris to compute the specific location of satellite
*/



/*
����ʱ��t
d240101//2024��1��1�գ�������һ
tk �Ӳο���Ԫ��ʼ���ʱ��, tk = t - toe, tΪ�۲�ʱ��
tyear
tmonth
tday
thour
tsec
*/
double toc_compute(double toe,double TOC_H, double TOC_Min, double TOC_Sec)
{
    double d240219 = 1;//��Ϊ������һ������1
    //double tyear = TOC_Y * 86400 * 30.6001 * 362.25;
   // double tmonth = TOC_M * 86400 * 30.6001;
  //  double tday = TOC_D * 86400;
    double thour = TOC_H * 3600;
    double tmin = TOC_Min * 60;
    double tsec = TOC_Sec;
    //toc�������������Ӳ�Ĳο�ʱ��
    double toc = d240219 * 86400 + thour + tmin + tsec;
    return toc;
}


/*
This is for computing the average angular velocity of satellite.
Input:
������
double sqrtA;//������ƽ����
������
double TOE;//�����Ĳο�ʱ�̣�GPS�����룩
��2��
double delta_n;//����ƽ���˶����������ֵ�Ĳ�ֵ����λ��rad/s����n�����˹�������صĳ����㶯
*/
double computeAverageAngularVelocity(double squareRootA, double delta_n)
{
	//�����˶���ƽ�����ٶ�
	double n = 0;
	//�ο�ʱ��Toe��ƽ�����ٶ�
	double n0 = 0;
	n0 = std::sqrt(GM) / (std::pow(squareRootA, 3));
	n = n0 + delta_n;
	return n;
}


/*
* This function is for computing the mean anomaly of satellite at the moment of observation.
* 	double M0;//�ο�ʱ���ƽ�����
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
����Բ�Ĳ�������x=acos�� �� y=bsin���У������Ǧȼ�Ϊƫ����ǡ�
Using rad to compute
input:
e and M
ʹ�õ����ⷨ��E�ĳ�ֵ������M����������1e-12������������1000
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
double omega;//���ص�Ǿ�,������
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
 ���ݿ���������ݣ�
 double IDOT;//i�ı仯��
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
    	double TOE;//�����Ĳο�ʱ�̣�GPS�����룩
    double deltaomega;//������ྭ�仯��(rad)
        double OMEGA;
        //ע�⣺�㲥�����и����Ĳ����ǲο�ʱ�̵�toeʱ��������ྭ�����Ǹ�ֵ�뱾����ʼʱ�̵ĸ������κ���ʱGAST_week֮��
        Ҳ����
        OMEGA = OMEGA_toe- GAST_week

        t������Ϊ����ʱ��
*/
double computeMomentOMEGA(double tk, double deltaomega, double OMEGA, double TOE)
{
    //����۲�˲���������ྭ
    double L = OMEGA + (deltaomega - OMEGA_E) * tk - OMEGA_E * TOE;
    return L;
}


/*
����������˲ʱ��������ϵ�е�λ��
���룺������Ĵ�ؾ���L�����ƽ������i�������ڹ��ƽ���λ��x,y
�����X Y Z
*/
void solveSatillitePosition(double L, double i,double x, double y,double *X, double *Y, double *Z)
{
    *X = (x * std::cos(L) - y * std::cos(i) * std::sin(L));
    *Y =(x * std::sin(L) + y * std::cos(i) * std::cos(L));
    *Z = (y * std::sin(i));
}



/*
�ɹ㲥�������������Ӳ�
*/
double clockBiasCompute(nav_body nav_b, double tk)
{
    //�Ӳ����

    //1.��Ը�����tr��GPS���Ƿ�Բ�ι��������������ЧӦ�������
    //1.1����Eccentric Anomaly
    double n = computeAverageAngularVelocity(nav_b.sqrtA, nav_b.delta_n);
    double M = computeInstantMeanAnomaly(nav_b.M0, n, tk);
    double EA = computeEccentricAnomaly(M, nav_b.e);

    double r_c = F * nav_b.e * nav_b.sqrtA * sin(EA);

    //1.2����SV clock correction
    double t_sv = nav_b.sa0 + nav_b.sa1 * tk + nav_b.sa2 * tk * tk + r_c;

    //1,3�õ������Ӳ�
    double satclock = t_sv - nav_b.TGD;

    return satclock;
}


/*
���������Ķ�����
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
    //�������������
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

    //���е���������ÿ�ε������µ���
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


//ecef(xyz) To ����������ϵENU
//���룺��վ��ecef������xr,yr,zr, ������ecef������xs,ys,zs
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

//�����Ǹ߶Ƚ�
//A����λ��
//rΪ�����򾶣�AΪ���Ƿ�λ�ǣ�hΪ���ǵĸ߶Ƚ�
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


/*��������*/
// �����ת��
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

//����ĳ˷�
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

// ����������������λ����
std::vector<std::vector<double>> createIdentityMatrix(int n)
{
    std::vector<std::vector<double>> identity(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        identity[i][i] = 1.0;
    }
    return identity;
}

// ��˹-Լ����Ԫ������
std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix)
{
    int n = matrix.size();
    std::vector<std::vector<double>> result = createIdentityMatrix(n);

    // ����������������޸�ԭʼ����
    std::vector<std::vector<double>> A = matrix;

    // ��˹-Լ����Ԫ��
    for (int i = 0; i < n; ++i) {
        // ���Խ���Ԫ�ع�һ
        double divisor = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= divisor;
            result[i][j] /= divisor;
        }

        // �����˶Խ���Ԫ�����������Ԫ������
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
����������Э���������ϵ�е�λ��
*/

