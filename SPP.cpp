#include"myhead.h"
using namespace std;

template<typename T>
// Overloaded elevcut function for 1D input vectors
std::vector<double> elevcut(const std::vector<T>& input, const std::vector<double>& elevation, double elevation_cutoff = 15.0) {
	size_t size = input.size();

	// Create mask
	std::vector<double> mask(size);
	for (size_t i = 0; i < size; ++i) {
		mask[i] = std::abs(elevation[i]) < elevation_cutoff ? std::numeric_limits<double>::quiet_NaN() : 1.0;
	}

	// Apply mask to input vector
	std::vector<double> inputs(size, 0.0);
	for (size_t i = 0; i < size; ++i) {
		if (!std::isnan(mask[i])) {
			inputs[i] = input[i];
		}
		else {
			inputs[i] = std::numeric_limits<double>::quiet_NaN();
		}
	}

	// Clear NaN values and create the output vector
	std::vector<double> output;
	for (size_t i = 0; i < size; ++i) {
		if (!std::isnan(inputs[i])) {
			output.push_back(inputs[i]);
		}
	}

	return output;
}
// Overloaded elevcut function for NavData input vectors
std::vector<NavData> elevcut(const std::vector<NavData>& input, const std::vector<double>& elevation, double elevation_cutoff = 15.0)
{
	size_t size = input.size();

	// Create mask
	std::vector<double> mask(size);
	for (size_t i = 0; i < size; ++i) {
		mask[i] = std::abs(elevation[i]) < elevation_cutoff ? std::numeric_limits<double>::quiet_NaN() : 1.0;
	}

	// Apply mask to input vector
	std::vector<NavData> inputs;
	inputs.resize(size);
	for (size_t i = 0; i < size; ++i) {
		if (!std::isnan(mask[i])) {
			inputs[i].X = input[i].X;
			inputs[i].Y = input[i].Y;
			inputs[i].Z = input[i].Z;
			inputs[i].clockBias = input[i].clockBias;
		}
		else
		{
			inputs[i].X = std::numeric_limits<double>::quiet_NaN();
			inputs[i].Y = std::numeric_limits<double>::quiet_NaN();
			inputs[i].Z = std::numeric_limits<double>::quiet_NaN();
			inputs[i].clockBias = std::numeric_limits<double>::quiet_NaN();
		}
	}

	// Clear NaN values and create the output vector
	std::vector<NavData> output;
	for (size_t i = 0; i < size; ++i) {
		if ((!std::isnan(inputs[i].X)) &&
			(!std::isnan(inputs[i].Y)) &&
			(!std::isnan(inputs[i].Z)) &&
			(!std::isnan(inputs[i].clockBias)))
		{
			output.push_back(inputs[i]);
		}
	}

	return output;
}

// Calculate GPS position
// Inputs:
// SOD = Second of day
// PRN = PRN number
//nav_body = 所有的N文件
// index = navigation index
// ps = pseudorange
// Outputs:
//NavData 包含x y z 和clock bias
NavData satpos_xyz_sbias(nav_body nav_b, ODate odate,int PRN, double SOD, double ps)
{
	//需要求的时刻
	//2024 02 19 00 00（存在了Odata里面的Odate）里算出周内秒之后，再加上SOD，
	double t = 86400 * 1 + SOD;
	//根据这个要求的时刻，在g_vector里面找到对应的星历
	double tk = 0;
	//计算t
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
	n = computeAverageAngularVelocity(nav_b.sqrtA, nav_b.delta_n);
	M = computeInstantMeanAnomaly(nav_b.M0, n, tk);
	E = computeEccentricAnomaly(M, nav_b.e);
	f = solveTrueAnomaly(E, nav_b.e);
	u_pie = solveArgumentOfLatitude(nav_b.omega, f);
	solvePerturbationCorrection(nav_b.Cuc, nav_b.Cus, nav_b.Crc, nav_b.Crs, nav_b.Cic, nav_b.Cis, &delta_u, &delta_r, &delta_i, u_pie);
	perturbationCorrection(nav_b.sqrtA, &u, &r, &i, nav_b.IDOT, E, tk, nav_b.e, u_pie, nav_b.i0, delta_u, delta_r, delta_i);
	computexyInOrbitPlane(r, u, &x, &y);
	L = computeMomentOMEGA(tk, nav_b.deltaomega, nav_b.OMEGA, nav_b.TOE);
	solveSatillitePosition(L, i, x, y, &X, &Y, &Z);
	/*
	std::cout << "年月日时分秒" << nav_b.TOC_Y << " " << nav_b.TOC_M << " " << nav_b.TOC_D << " " << nav_b.TOC_H << " " << nav_b.TOC_Min << " " << nav_b.TOC_Sec << endl;
	std::cout << "卫星标识:" << nav_b.monument << std::endl;
	std::cout << "PRN是" << nav_b.sPRN << endl;
	std::cout << "计算得到的卫星坐标是" << endl;
	std::cout << "X:" << std::setprecision(9) << X << endl;
	std::cout << "Y:" << std::setprecision(9) << Y << endl;
	std::cout << "Z:" << std::setprecision(9) << Z << endl;
	*/
	NavData navData;
	navData.monument = nav_b.monument;
	navData.sPRN = nav_b.sPRN;
	navData.tk = tk;
	navData.X = X;
	navData.Y = Y;
	navData.Z = Z;
	return navData;
}


/*
计算欧氏距离
*/
// Function to compute Euclidean distance
std::vector<double> computeRefDistance(const rcvposCoor& refpos, const std::vector<NavData>& satpos_m)
{
	size_t numSatellites = satpos_m.size();
	std::vector<double> ctest(numSatellites);
	for (size_t i = 0; i < numSatellites; ++i) {
		double sum = 0.0;
		sum = std::pow(refpos.x - satpos_m[i].X, 2) + std::pow(refpos.y - satpos_m[i].Y, 2) + std::pow(refpos.z - satpos_m[i].Z, 2);
		ctest[i] = std::sqrt(sum);
	}
	return ctest;
}

// Function to remove outliers
void removeOutliers(std::vector<double>& cP1, std::vector<NavData>& satpos_m, const rcvposCoor& refpos)
{
	std::vector<double> ctest = computeRefDistance(refpos, satpos_m);
	std::vector<double> cError(cP1.size());
	std::transform(cP1.begin(), cP1.end(), ctest.begin(), cError.begin(), [](double p1, double test) {
		return p1 - test;
		});
	// Remove elements where the absolute error is greater than 1000
	std::vector<double>::iterator itP1 = cP1.begin();
	std::vector<NavData>::iterator itSat = satpos_m.begin();
	for (auto itErr = cError.begin(); itErr != cError.end();) {
		if (std::abs(*itErr) > 1000) {
			itP1 = cP1.erase(itP1);
			itSat = satpos_m.erase(itSat);
			itErr = cError.erase(itErr);
		}
		else {
			++itP1;
			++itSat;
			++itErr;
		}
	}

}


vector<NavData> FlightTimeCorrection(vector<NavData>satpos,vector<double> dTflightSeconds)
{

	vector<NavData> newsatpos;
	newsatpos.resize(satpos.size());
	for (int i = 0;i<satpos.size();i++)
	{
		double theta = OMEGA_E * dTflightSeconds[i];
		newsatpos[i].X = cos(theta) * satpos[i].X + sin(theta) * satpos[i].Y;
		newsatpos[i].Y = -sin(theta) * satpos[i].X + cos(theta) * satpos[i].Y;
		newsatpos[i].Z = satpos[i].Z;
	}
	return newsatpos;
}




//开始单点定位
//输入：
//1.obs：o文件读好的数据
//2.pnavHead: N文件中头文件的参数，
//3.g_vector：一个向量，里面包含了N文件中所有GPS卫星的参考星历
//4.Spath：
void positioning(const OData& obs, vector<nav_body> g_vector,pnav_head pnavHead)
{

	//设置好需要导出数据所在的文件名
	string filename = "SPP_result.csv";
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "无法打开文件: " << filename << std::endl;
		return;
	}
	// 写入CSV表头
	file << "SOD,Num_Sat_Rcv, X, Y, Z, Receiver_Bias,GDOP" << std::endl;
	//设置好相关参数
	double stt = 0;//起始的UTC时间，单位为hour
	double stp = 25;//结束的UTC时间，单位为hour
	int interval= 30; //每次单点定位间隔的秒数，因为O文件30秒新增一项

	//矩阵长度
	int len_matrix = (int)((stp - stt) * 3600 - 1);//比如25*3600-1 = 8999

	//将rcv.pos输入作为refpos
	ECEF refpos;
	refpos.x = obs.rcvpos[0].x;
	refpos.y = obs.rcvpos[0].y;
	refpos.z = obs.rcvpos[0].z;
	//将refpos（地心地固）转化为BLH
	BLH refpos_blh = xyz2blh(refpos);
	std::cout << "转化后的refpos是" << refpos_blh.B<<" " << refpos_blh.L << " " << refpos_blh.H << endl;

	//初始化单点定位结果的结构体
	//这里为什么是32？其实是因为GPS有32颗卫星
	// 初始化位置结果结构体
	PositioningResult userpos;
	userpos.blh.resize(len_matrix, { 0, 0, 0 });
	userpos.xyz.resize(len_matrix, { 0, 0, 0 });
	userpos.bs.resize(len_matrix, std::vector<double>(32, 0));
	userpos.br.resize(len_matrix, 0);
	userpos.elevation.resize(len_matrix, std::vector<double>(32, 0));

	//初始化举例误差结构体
	DistError disterr;
	disterr.horizontal.resize(len_matrix, NAN);
	disterr.EW.resize(len_matrix, NAN);
	disterr.NS.resize(len_matrix, NAN);
	disterr.height.resize(len_matrix, NAN);

	//初始化模型结构体
	Model model;
	model.tropo.resize(len_matrix, std::vector<double>(32, NAN));
	model.iono.resize(len_matrix, std::vector<double>(32, NAN));
	model.hdop.resize(len_matrix, NAN);

	//第一步，初始化矩阵G
	std::vector<std::vector<double>> G(4, std::vector<double>(4, 0));

	//GPS卫星单颗卫星数据，分32个存储
	vector<vector<nav_body>>  gpsSats;
	gpsSats.resize(32);
	//Read Body File.

	std::cout << "完毕" << endl;
	//将时间转化为gps时
	std::vector<NavData> data;

	for (nav_body& nav_b : g_vector)
	{
		gpsSats[nav_b.sPRN - 1].push_back(nav_b);
	}

	//开始计算
	for (int mode = 4; mode <= 4; ++mode)
	{
		//Initialize Uer position
		Vector4D x0;
		x0.x = 0; x0.y = 0; x0.z =0; x0.br = 0;
		Vector3D xyz0 = {x0.x,x0.y,x0.z};//提取向量 x0 的前 3 个元素
		//这里x0,y0,z0是迭代的初值，是用户位置的近似值
		double br = x0.br;

		for (double i = (stt * 3600) + 1; i < (stp * 3600) - 1; i += interval)
		{

			//找到同一SOD内对应的所有符合条件的index，存储在ephIndexes里
			vector<int> ephIndexes;
			auto it = find(obs.epoch.begin(), obs.epoch.end(), i - 1);
			while (it != obs.epoch.end())
			{
				int index = distance(obs.epoch.begin(), it);
				ephIndexes.push_back(index);
				it = find(it + 1, obs.epoch.end(), i - 1);
			}

			//根据储存的索引，按索引提取出卫星的编号，储存在PRNs里，卫星编号是可以重复的
			vector<int>PRNs;
			for (int index : ephIndexes)
			{
				PRNs.push_back(obs.index[index]);
			}

			//如果这个时间点观测到的卫星数少于四个，则跳过
			if (PRNs.size() < 4)
				continue;


			//Second Of Day, SOD的计算公式为：日内秒+（小时数*3600+分钟数*60+自带秒数）
			//括号里面是参考时刻的秒数

			vector<double> SOD;
			for (int index : ephIndexes)
			{
				double sod = obs.epoch[index] + obs.date.hour * 3600 + obs.date.minute * 60 + obs.date.second;
				SOD.push_back(sod);
			}

			//第二步，从O文件读取伪距观测值
			//Pseudorange
			vector<double>C1;
			for (int index : ephIndexes)
			{
				double c1 = obs.data[index][0];
				C1.push_back(c1);
			}

			//第三步，计算卫星位置和卫星钟差
			vector<NavData> satpos_vector;
			satpos_vector.resize(PRNs.size());
			//按照PRNs里的编号把卫星位置存起来
			//-----------------------------这一步等到卫星位置计算的内容进一步完成之后，再补上
			for (int index = 0;index<PRNs.size();index++)
			{
				double t = 86400 + SOD[index];
				//这个t是信号到接收机的时间，卫星发射信号的时间还得减去一个tr = 伪距/光速.
				double tr = C1[index] / C;
				nav_body trueParams = findTrueOne(t, gpsSats[PRNs[index]-1]);
				double tk = t - tr - trueParams.TOE;
				//计算卫星的位置以及此时卫星的钟差
			    NavData navdata = myComputeLocationOfSatellite(trueParams, tk);
				satpos_vector[index] = navdata;
			}
			//Check number of satellite
			int OK_num = 0;
			for (NavData satpos : satpos_vector)
			{
				if (satpos.monument == "G")
					OK_num++;
			}
			if (OK_num < 4)
				continue;
			//第四步，计算用户位置
			//首先重新设置好改正模型矩阵
			vector<double> dIon_klob;
			dIon_klob.resize(PRNs.size(), NAN);
			vector<double>tropoDelay;
			tropoDelay.resize(PRNs.size(), NAN);
			vector<double>ionoDelay;
			ionoDelay.resize(PRNs.size(), NAN);

			int stop = 0;
			double differenceOfPositioning = 10e-4;
			//定义迭代次数
			int iterationTime = 0;
			double dxyz = 99999;
			while (dxyz > differenceOfPositioning)
			{
				stop++;
				//最多迭代29次
				if (stop == 30)
					break;
				vector<RAH> rah_vector;
				ECEF ecefxyz0;
				ecefxyz0.x = xyz0.x; ecefxyz0.y = xyz0.y; ecefxyz0.z = xyz0.z;
				BLH userpos_blh = xyz2blh(ecefxyz0);
				ECEF receiver;
				receiver.x = xyz0.x; receiver.y = xyz0.y; receiver.z = xyz0.z;
				for (int i = 0; i < satpos_vector.size(); i++)
				{
					ECEF satellite = { satpos_vector[i].X,satpos_vector[i].Y,satpos_vector[i].Z };
					RAH rah;
					rah = sateliite_rah(receiver, satellite);
					rah_vector.push_back(rah);
					//计算对流层延迟
					double d_hyd = 0;
					double d_wet = 0;
					EstimateTropDalay(userpos_blh.L, userpos_blh.H, 50, d_hyd, d_wet);
					double delaytropo = (d_hyd + d_wet) * (1.001 / sqrt(0.002001 + sin(rah.H*PI/180) * sin(rah.H * PI / 180.0)  ));
					tropoDelay[i] = delaytropo;
				}

				double ionprm[2][4] = {
				 {pnavHead->ION_alpha[0],pnavHead->ION_alpha[1],pnavHead->ION_alpha[2],pnavHead->ION_alpha[3]},
				 {pnavHead->ION_beta[0],pnavHead->ION_beta[1],pnavHead->ION_beta[2],pnavHead->ION_beta[3]}
				};

				for (int ii = 0; ii < PRNs.size(); ii++)
				{
					dIon_klob[ii] = klobuchar_model(userpos_blh.L, userpos_blh.B, rah_vector[ii].H, SOD[ii], ionprm);
					ionoDelay[ii] = C * dIon_klob[ii];
				}

				//把所有高度角存在一个向量里
				std::vector<double> elevs;
				elevs.resize(rah_vector.size()); // 预留足够的空间以避免重复分配
				for (int i = 0; i <elevs.size(); i++)
					elevs[i] = rah_vector[i].H;
				//高度角 cutoff 15 degrees，因为高度角小于15的，都不准确
				double elev_mask = 15.0;

				std::vector <double>elev_m;
				std::vector<double>C1_m;
				std::vector<NavData>satpos_m;
				std::vector<double>d_tropo;
				std::vector<double>d_iono;
				std::vector<double>PRN_m;
				if (stop > 1)
				{
					elev_m = elevcut(elevs, elevs, elev_mask);
					C1_m = elevcut(C1, elevs, elev_mask);
					satpos_m = elevcut(satpos_vector, elevs, elev_mask);
					d_tropo = elevcut(tropoDelay, elevs, elev_mask);
					d_iono = elevcut(ionoDelay, elevs, elev_mask);
					PRN_m = elevcut(PRNs, elevs, elev_mask);
				}
				else
				{
					std::transform(elevs.begin(), elevs.end(), tropoDelay.begin(), [](double value) {return value * 0.0; });
					std::transform(elevs.begin(), elevs.end(), ionoDelay.begin(), [](double value) {return value * 0.0; });
					elev_m = elevcut(elevs, elevs, 0);
					C1_m = elevcut(C1, elevs, 0);
					satpos_m = elevcut(satpos_vector, elevs, 0);
					d_tropo = elevcut(tropoDelay, elevs, 0);
					d_iono = elevcut(ionoDelay, elevs, 0);
					PRN_m = elevcut(PRNs, elevs, 0);
				}
				vector<double>cP1;
				cP1.resize(C1_m.size(), NAN);

				//Perform Pseudorange Correction
				if (cP1.size() != satpos_m.size())
					break;//如果伪距数量不等于卫星的数量，则该部分数据无法用于计算

				switch (mode)
				{
					case NO_CORRECTION:
						for (size_t i = 0; i < cP1.size(); ++i) {
							cP1[i] = C1_m[i] + C * satpos_m[i].clockBias;
						}
						break;
					case IONO_CORRECTION:
						for (size_t i = 0; i < cP1.size(); ++i) {
							cP1[i] = C1_m[i] + C * satpos_m[i].clockBias - d_iono[i];
						}
						break;
					case TROPO_CORRECTION:
						for (size_t i = 0; i < cP1.size(); ++i) {
							cP1[i] = C1_m[i] + C * satpos_m[i].clockBias - d_tropo[i];
						}
						break;
					case TROPO_IONO_CORRECTION:
						for (size_t i = 0; i < cP1.size(); ++i) {
							cP1[i] = C1_m[i] + C * satpos_m[i].clockBias - d_tropo[i] - d_iono[i];
						}
						break;
					default:
						std::cerr << "Invalid mode!" << std::endl;
						break;
				}

				//remove outlier，移除掉与真实距离相差超过1000米的伪距
				//removeOutliers(cP1, satpos_m, obs.rcvpos[0]);

				if (satpos_m.size() < 4 || C1_m.size() < 4 || satpos_m.size() != C1_m.size())
					continue;


				vector<double>dtflights(satpos_m.size());
				//用br和clock error对卫星航向进行改正
				for (int i = 0; i < satpos_m.size(); i++)
				{
					double dtflight = (cP1[i] - br) / C + satpos_m[i].clockBias;
					dtflights[i] = dtflight;
				}
				vector<NavData>satpos_cor(satpos_m.size());
				satpos_cor = FlightTimeCorrection(satpos_m, dtflights);

				//计算卫星到初值x0的向量和距离
				vector<Vector3D> v;
				v.resize(satpos_cor.size());
				for (int i = 0; i < v.size(); i++)
				{
					v[i].x = xyz0.x - satpos_cor[i].X;
					v[i].y = xyz0.y - satpos_cor[i].Y;
					v[i].z = xyz0.z - satpos_cor[i].Z;
				}
				vector<double> range_v;
				double sum = 0;
				for (int i = 0; i < v.size(); i++)
				{
					sum = v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z;
					double range = sqrt(sum);
					range_v.push_back(range);
				}

				
				for (int i = 0; i < v.size(); i++)
				{
					v[i].x = v[i].x / range_v[i];
					v[i].y = v[i].y / range_v[i];
					v[i].z = v[i].z / range_v[i];
				}

				//计算a-priori range residual
				vector<double> prHat(range_v.size());
				std::transform(range_v.begin(), range_v.end(), prHat.begin(), [br](double i) {return i + br; });
				vector<double>P(prHat.size());
				for (int i = 0; i < P.size(); i++)
					P[i] = cP1[i] - prHat[i]; //P阵

			    //创建H阵，H阵的第一列是V.T，第二列是1
				//也就是x1 y1 z1 1
				//             x2 y2 z2 1
				//             ...
				//             xm ym zm 1
				vector<vector<double>>H(P.size(),std::vector<double>(4,0.0));
				for (int i = 0; i < H.size();i++)
				{
					H[i][0] = v[i].x;
					H[i][1] = v[i].y;
					H[i][2] = v[i].z;
					H[i][3] = 1;
				}

				vector<vector<double>>P_matrix(P.size(), std::vector<double>(1, 0.0));
				for (int i = 0; i < P_matrix.size(); i++)
					P_matrix[i][0] = P[i];

					


				//LSE solution
				// 计算H的转置矩阵
				std::vector<std::vector<double>> H_transpose = transpose(H);
				// 计算H'*H      4*m    x    m*4
				std::vector<std::vector<double>> H_transpose_multiply_H = matrixMultiply(H_transpose, H);
				// 计算H'*P      4*m    x    m*1
				std::vector<std::vector<double>> H_transpose_multiply_P = matrixMultiply(H_transpose, P_matrix);
				// 计算(H'*H)的逆
				std::vector<std::vector<double>> H_transpose_multiply_H_inverse = inverse(H_transpose_multiply_H);
				// 计算dx = inv(H'*H)*H'*P
				std::vector<std::vector<double>>dx = matrixMultiply(H_transpose_multiply_H_inverse, H_transpose_multiply_P);

				//更新xyz0和br
				xyz0.x += dx[0][0];
				xyz0.y += dx[1][0];
				xyz0.z += dx[2][0];
				br += dx[3][0];
				G = H_transpose_multiply_H_inverse;
				dxyz = pow(dx[0][0], 2) + pow(dx[1][0], 2) + pow(dx[2][0], 2) + pow(dx[3][0], 2);
				iterationTime++;
				std::cout << "已经迭代了" << iterationTime << "次" << endl;
			}
			double GDOP = sqrt(G[0][0] + G[1][1]+ G[2][2] + G[3][3]);
			std::cout << "现在是2月19日，日内秒" << i << "秒" << endl;
			std::cout<<"收到" << PRNs.size()<<"颗GPS卫星的伪距信号，他们是：" << endl;
			for (int prn : PRNs)
				std::cout <<"G" << prn << " ";
			std::cout << endl << "这是第" << mode << "种模式，单点定位得到的用户的坐标为" << endl;
			std::cout<<std::setprecision(9) << "x是" << xyz0.x << ",y是" << xyz0.y << ",z是" << xyz0.z << ",接收机钟差是" << br << endl;
			std::cout << "-----------------------------------------分割线----------------------------------------" << endl;
			//	file << "SOD,Num_Sat_Rcv, X, Y, Z, Receiver_Bias,GDOP" << std::endl;
			file << i << ","
				<< PRNs.size()<< ","
				<< std::setprecision(9) << xyz0.x << ","
				<< std::setprecision(9) << xyz0.y << ","
				<< std::setprecision(9) << xyz0.z << ","
				<< std::setprecision(9) << br << ","
				<<std::setprecision(9)<<GDOP<<std::endl;
		}
	}
	file.close();
	std::cout << "单点定位结果已经保存到了文件"<<filename << std::endl;
}


//开始无电离层组合单点定位
//输入：
//1.obs：o文件读好的数据
//2.pnavHead: N文件中头文件的参数，
//3.g_vector：一个向量，里面包含了N文件中所有GPS卫星的参考星历
void positioning_iono_free(const OData& obs, vector<nav_body> g_vector, pnav_head pnavHead)
{
	//设置好需要导出数据所在的文件名
	string filename = "SPP_result_using_IF_combination.csv";
	std::ofstream file(filename);
	if (!file.is_open()) 
	{
		std::cerr << "无法打开文件: " << filename << std::endl;
		return;
	}
	// 写入CSV表头
	file << "SOD,Num_Sat_Rcv, X, Y, Z, Receiver_Bias,GDOP" << std::endl;
	//设置好相关参数
	double stt = 0;//起始的UTC时间，单位为hour
	double stp = 25;//结束的UTC时间，单位为hour
	int interval = 30; //每次单点定位间隔的秒数，因为O文件30秒新增一项

	//矩阵长度
	int len_matrix = (int)((stp - stt) * 3600 - 1);//比如25*3600-1 = 8999

	//将rcv.pos输入作为refpos
	ECEF refpos;
	refpos.x = obs.rcvpos[0].x;
	refpos.y = obs.rcvpos[0].y;
	refpos.z = obs.rcvpos[0].z;
	//将refpos（地心地固）转化为BLH
	BLH refpos_blh = xyz2blh(refpos);
	std::cout << "转化后的refpos是" << refpos_blh.B << " " << refpos_blh.L << " " << refpos_blh.H << endl;

	//初始化单点定位结果的结构体
	//这里为什么是32？其实是因为GPS有32颗卫星
	// 初始化位置结果结构体
	PositioningResult userpos;
	userpos.blh.resize(len_matrix, { 0, 0, 0 });
	userpos.xyz.resize(len_matrix, { 0, 0, 0 });
	userpos.bs.resize(len_matrix, std::vector<double>(32, 0));
	userpos.br.resize(len_matrix, 0);
	userpos.elevation.resize(len_matrix, std::vector<double>(32, 0));

	//初始化举例误差结构体
	DistError disterr;
	disterr.horizontal.resize(len_matrix, NAN);
	disterr.EW.resize(len_matrix, NAN);
	disterr.NS.resize(len_matrix, NAN);
	disterr.height.resize(len_matrix, NAN);

	//初始化模型结构体
	Model model;
	model.tropo.resize(len_matrix, std::vector<double>(32, NAN));
	model.hdop.resize(len_matrix, NAN);

	//第一步，初始化矩阵G
	std::vector<std::vector<double>> G(4, std::vector<double>(4, 0));

	//GPS卫星单颗卫星数据，分32个存储
	vector<vector<nav_body>>  gpsSats;
	gpsSats.resize(32);
	//Read Body File.

	std::cout << "完毕" << endl;
	//将时间转化为gps时
	std::vector<NavData> data;

	for (nav_body& nav_b : g_vector)
	{
		gpsSats[nav_b.sPRN - 1].push_back(nav_b);
	}

	//开始计算
			//Initialize Uer position
	Vector4D x0;
	x0.x = 0; x0.y = 0; x0.z = 0; x0.br = 0;
	Vector3D xyz0 = { x0.x,x0.y,x0.z };//提取向量 x0 的前 3 个元素
	//这里x0,y0,z0是迭代的初值，是用户位置的近似值
	double br = x0.br;

	for (double i = (stt * 3600) + 1; i < (stp * 3600) - 1; i += interval)
	{

		//找到同一SOD内对应的所有符合条件的index，存储在ephIndexes里
		vector<int> ephIndexes;
		auto it = find(obs.epoch.begin(), obs.epoch.end(), i - 1);
		while (it != obs.epoch.end())
		{
			int index = distance(obs.epoch.begin(), it);
			ephIndexes.push_back(index);
			it = find(it + 1, obs.epoch.end(), i - 1);
		}
		//根据储存的索引，按索引提取出卫星的编号，储存在PRNs里，卫星编号是可以重复的
		vector<int>PRNs;
		for (int index : ephIndexes)
		{
			PRNs.push_back(obs.index[index]);
		}

		//如果这个时间点观测到的卫星数少于四个，则跳过
		if (PRNs.size() < 4)
			continue;


		//Second Of Day, SOD的计算公式为：日内秒+（小时数*3600+分钟数*60+自带秒数）
		//括号里面是参考时刻的秒数

		vector<double> SOD;
		for (int index : ephIndexes)
		{
			double sod = obs.epoch[index] + obs.date.hour * 3600 + obs.date.minute * 60 + obs.date.second;
			SOD.push_back(sod);
		}

		//第二步，从O文件读取伪距观测值
		//读取C1C观测值
		vector<double>C1C;
		//C2X Values are stored in this vector.
		vector<double>C2X;
		for (int index : ephIndexes)
		{
			C1C.push_back(obs.data[index][0]);
			C2X.push_back(obs.data[index][1]);
		}

		//Calculate the linear Combinations of Observations.
		vector<double>IF_mn;
		for (int i = 0;i<C1C.size();i++)
		{
			double if_mn = 2.54573 * C1C[i] - 1.54573 * C2X[i];
			IF_mn.push_back(if_mn);
		}


		//第三步，计算卫星位置和卫星钟差
		vector<NavData> satpos_vector;
		satpos_vector.resize(PRNs.size());
		//按照PRNs里的编号把卫星位置存起来
		//-----------------------------这一步等到卫星位置计算的内容进一步完成之后，再补上
		for (int index = 0; index < PRNs.size(); index++)
		{
			double t = 86400 + SOD[index];
			//这个t是信号到接收机的时间，卫星发射信号的时间还得减去一个tr = 伪距/光速.
			double tr = IF_mn[index] / C;
			nav_body trueParams = findTrueOne(t, gpsSats[PRNs[index] - 1]);
			double tk = t - tr - trueParams.TOE;
			//计算卫星的位置以及此时卫星的钟差
			NavData navdata = myComputeLocationOfSatellite(trueParams, tk);
			satpos_vector[index] = navdata;
		}
		//Check number of satellite
		int OK_num = 0;
		for (NavData satpos : satpos_vector)
		{
			if (satpos.monument == "G")
				OK_num++;
		}
		if (OK_num < 4)
			continue;
		//第四步，计算用户位置
		//首先重新设置好改正模型矩阵
		vector<double> dIon_klob;
		dIon_klob.resize(PRNs.size(), NAN);
		vector<double>tropoDelay;
		tropoDelay.resize(PRNs.size(), NAN);
		vector<double>ionoDelay;
		ionoDelay.resize(PRNs.size(), NAN);

		int stop = 0;
		double differenceOfPositioning = 10e-4;
		//定义迭代次数
		int iterationTime = 0;
		double dxyz = 99999;
		while (dxyz > differenceOfPositioning)
		{
			stop++;
			//最多迭代29次
			if (stop == 30)
				break;
			vector<RAH> rah_vector;
			ECEF ecefxyz0;
			ecefxyz0.x = xyz0.x; ecefxyz0.y = xyz0.y; ecefxyz0.z = xyz0.z;
			BLH userpos_blh = xyz2blh(ecefxyz0);
			ECEF receiver;
			receiver.x = xyz0.x; receiver.y = xyz0.y; receiver.z = xyz0.z;
			for (int i = 0; i < satpos_vector.size(); i++)
			{
				ECEF satellite = { satpos_vector[i].X,satpos_vector[i].Y,satpos_vector[i].Z };
				RAH rah;
				rah = sateliite_rah(receiver, satellite);
				rah_vector.push_back(rah);
				//计算对流层延迟
				double d_hyd = 0;
				double d_wet = 0;
				EstimateTropDalay(userpos_blh.L, userpos_blh.H, 50, d_hyd, d_wet);
				double delaytropo = (d_hyd + d_wet) * (1.001 / sqrt(0.002001 + sin(rah.H * PI / 180) * sin(rah.H * PI / 180.0)));
				tropoDelay[i] = delaytropo;
			}

			//把所有高度角存在一个向量里
			std::vector<double> elevs;
			elevs.resize(rah_vector.size()); // 预留足够的空间以避免重复分配
			for (int i = 0; i < elevs.size(); i++)
				elevs[i] = rah_vector[i].H;
			//高度角 cutoff 15 degrees，因为高度角小于15的，都不准确
			double elev_mask = 15.0;

			std::vector <double>elev_m;
			std::vector<double>IF_mn_m;
			std::vector<NavData>satpos_m;
			std::vector<double>d_tropo;
			std::vector<double>PRN_m;
			if (stop > 1)
			{
				elev_m = elevcut(elevs, elevs, elev_mask);
				IF_mn_m = elevcut(IF_mn, elevs, elev_mask);
				satpos_m = elevcut(satpos_vector, elevs, elev_mask);
				d_tropo = elevcut(tropoDelay, elevs, elev_mask);
				PRN_m = elevcut(PRNs, elevs, elev_mask);
			}
			else
			{
				std::transform(elevs.begin(), elevs.end(), tropoDelay.begin(), [](double value) {return value * 0.0; });
				std::transform(elevs.begin(), elevs.end(), ionoDelay.begin(), [](double value) {return value * 0.0; });
				elev_m = elevcut(elevs, elevs, 0);
				IF_mn_m = elevcut(IF_mn, elevs, 0);
				satpos_m = elevcut(satpos_vector, elevs, 0);
				d_tropo = elevcut(tropoDelay, elevs, 0);
				PRN_m = elevcut(PRNs, elevs, 0);
			}
			vector<double>IF_mn_P1;
			IF_mn_P1.resize(IF_mn_m.size(), NAN);

			//Perform Pseudorange Correction
			if (IF_mn_P1.size() != satpos_m.size())
				break;//如果伪距数量不等于卫星的数量，则该部分数据无法用于计算
			for (size_t i = 0; i < IF_mn_P1.size(); ++i) 
				IF_mn_P1[i] = IF_mn_m[i] + C * satpos_m[i].clockBias - d_tropo[i];
				//remove outlier，移除掉与真实距离相差超过1000米的伪距
				//removeOutliers(cP1, satpos_m, obs.rcvpos[0]);

				if (satpos_m.size() < 4 || IF_mn_m.size() < 4 || satpos_m.size() != IF_mn_m.size())
					continue;


				vector<double>dtflights(satpos_m.size());
				//用br和clock error对卫星航向进行改正
				for (int i = 0; i < satpos_m.size(); i++)
				{
					double dtflight = (IF_mn_P1[i] - br) / C + satpos_m[i].clockBias;
					dtflights[i] = dtflight;
				}
				vector<NavData>satpos_cor(satpos_m.size());
				satpos_cor = FlightTimeCorrection(satpos_m, dtflights);

				//计算卫星到初值x0的向量和距离
				vector<Vector3D> v;
				v.resize(satpos_cor.size());
				for (int i = 0; i < v.size(); i++)
				{
					v[i].x = xyz0.x - satpos_cor[i].X;
					v[i].y = xyz0.y - satpos_cor[i].Y;
					v[i].z = xyz0.z - satpos_cor[i].Z;
				}
				vector<double> range_v;
				double sum = 0;
				for (int i = 0; i < v.size(); i++)
				{
					sum = v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z;
					double range = sqrt(sum);
					range_v.push_back(range);
				}


				for (int i = 0; i < v.size(); i++)
				{
					v[i].x = v[i].x / range_v[i];
					v[i].y = v[i].y / range_v[i];
					v[i].z = v[i].z / range_v[i];
				}

				//计算a-priori range residual
				vector<double> prHat(range_v.size());
				std::transform(range_v.begin(), range_v.end(), prHat.begin(), [br](double i) {return i + br; });
				vector<double>P(prHat.size());
				for (int i = 0; i < P.size(); i++)
					P[i] = IF_mn_P1[i] - prHat[i]; //P阵

				//创建H阵，H阵的第一列是V.T，第二列是1
				//也就是x1 y1 z1 1
				//             x2 y2 z2 1
				//             ...
				//             xm ym zm 1
				vector<vector<double>>H(P.size(), std::vector<double>(4, 0.0));
				for (int i = 0; i < H.size(); i++)
				{
					H[i][0] = v[i].x;
					H[i][1] = v[i].y;
					H[i][2] = v[i].z;
					H[i][3] = 1;
				}

				vector<vector<double>>P_matrix(P.size(), std::vector<double>(1, 0.0));
				for (int i = 0; i < P_matrix.size(); i++)
					P_matrix[i][0] = P[i];




				//LSE solution
				// 计算H的转置矩阵
				std::vector<std::vector<double>> H_transpose = transpose(H);
				// 计算H'*H      4*m    x    m*4
				std::vector<std::vector<double>> H_transpose_multiply_H = matrixMultiply(H_transpose, H);
				// 计算H'*P      4*m    x    m*1
				std::vector<std::vector<double>> H_transpose_multiply_P = matrixMultiply(H_transpose, P_matrix);
				// 计算(H'*H)的逆
				std::vector<std::vector<double>> H_transpose_multiply_H_inverse = inverse(H_transpose_multiply_H);
				// 计算dx = inv(H'*H)*H'*P
				std::vector<std::vector<double>>dx = matrixMultiply(H_transpose_multiply_H_inverse, H_transpose_multiply_P);

				//更新xyz0和br
				xyz0.x += dx[0][0];
				xyz0.y += dx[1][0];
				xyz0.z += dx[2][0];
				br += dx[3][0];
				G = H_transpose_multiply_H_inverse;
				dxyz = pow(dx[0][0], 2) + pow(dx[1][0], 2) + pow(dx[2][0], 2) + pow(dx[3][0], 2);
				iterationTime++;
				std::cout << "已经迭代了" << iterationTime << "次" << endl;
			}
			double GDOP = sqrt(G[0][0] + G[1][1] + G[2][2] + G[3][3]);
			std::cout << "现在是2月19日，日内秒" << i-1 << "秒" << endl;
			std::cout << "收到" << PRNs.size() << "颗GPS卫星的伪距信号，他们是：" << endl;
			for (int prn : PRNs)
				std::cout << "G" << prn << " ";
			std::cout << endl << "这是第无电离层组合模式，单点定位得到的用户的坐标为" << endl;
			std::cout << std::setprecision(11) << "x是" << xyz0.x << ",y是" << xyz0.y << ",z是" << xyz0.z << ",接收机钟差是" << br << endl;
			std::cout << std::setprecision(11) <<"GDOP是" << GDOP <<std::endl;
			std::cout << "-----------------------------------------分割线----------------------------------------" << endl;
			file << i << ","
				<< PRNs.size() << ","
				<< std::setprecision(11) << xyz0.x << ","
				<< std::setprecision(11) << xyz0.y << ","
				<< std::setprecision(11) << xyz0.z << ","
				<< std::setprecision(11) << br << ","
				<< std::setprecision(11) << GDOP << std::endl;
		}
	file.close();
	std::cout << "单点定位结果已经保存到了文件" << filename << std::endl;
}