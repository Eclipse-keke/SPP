#include"myhead.h"



//电离层改正模型
// Ionospheric delay calculation(Klobuchar model, 1987)
// Inputs:
// fi = latitude of rcv position
// lambda = longitude of rcv position
// elev = elevation angle
// tow = time of week
// ionprm = navigation coeficients(alpha and beta)
// Outputs:
// dIon = Ionospheric delay
double klobuchar_model(double fi, double lambda, double elev, double tow, const double ionprm[2][4]) 
{
	const double* alfa = ionprm[0];  // alpha
	const double* beta = ionprm[1];  // Beta

	double azimuth = 360;
	double c = 2.99792458e8;              // speed of light
	double deg2semi = 1.0 / 180.0;        // degrees to semicircles
	double semi2rad = PI;               // semicircles to radians
	double deg2rad = PI / 180.0;        // degrees to radians

	double a = azimuth * deg2rad;         // azimuth in radians
	double e = elev * deg2semi;           // elevation angle in semicircles

	double psi = 0.0137 / (e + 0.11) - 0.022; // Earth Centered angle

	double lat_i = fi * deg2semi + psi * cos(a);  // Subionospheric lat
	if (lat_i > 0.416)
		lat_i = 0.416;
	else if (lat_i < -0.416)
		lat_i = -0.416;

	// Subionospheric long
	double long_i = lambda * deg2semi + (psi * sin(a) / cos(lat_i * semi2rad));

	// Geomagnetic latitude
	double lat_m = lat_i + 0.064 * cos((long_i - 1.617) * semi2rad);

	double t = 4.32e4 * long_i + tow;
	t = fmod(t, 86400.0);  // Seconds of day
	if (t > 86400.0)
		t -= 86400.0;
	if (t < 0.0)
		t += 86400.0;

	double sF = 1.0 + 16.0 * pow((0.53 - e), 3);  // Slant factor

	// Period of model
	double PER = beta[0] + beta[1] * lat_m + beta[2] * pow(lat_m, 2) + beta[3] * pow(lat_m, 3);

	if (PER < 72000.0)
		PER = 72000.0;

	double x = 2.0 * PI * (t - 50400.0) / PER;  // Phase of the model

	// Amplitude of the model
	double AMP = alfa[0] + alfa[1] * lat_m + alfa[2] * pow(lat_m, 2) + alfa[3] * pow(lat_m, 3);
	if (AMP < 0.0)
		AMP = 0.0;

	// Ionospheric correction
	double dIon1;
	if (abs(x) > 1.57)
		dIon1 = sF * 5.0e-9;
	else
		dIon1 = sF * (5.0e-9 + AMP * (1.0 - x * x / 2.0 + pow(x, 4) / 24.0));
	return dIon1;  // * 10^9
}


// Function for linear interpolation
double Interpol(double Wp, double Wmax, double Wmin, double Vmax, double Vmin) {
	return Vmin + (Vmax - Vmin) * ((Wp - Wmin) / (Wmax - Wmin));
}

/*
对流层改正模型
*/
void EstimateTropDalay(double Latitude, double Height, double DOY, double& d_hyd, double& d_wet) 
{
	// Check for NaN values
	if (std::isnan(Latitude) || std::isnan(DOY)) {
		throw std::invalid_argument("Value of Latitude and DOY cannot be NaN");
	}

	double Lat = std::abs(Latitude);
	int Dmin = (Latitude >= 0) ? 28 : 211;

	if (Height >= 1000) Height = 1000; // clear outage
	double H = Height;

	// Constants
	double k1 = 77.604; // K/mbar
	double k2 = 382000; // K^2/mbar
	double Rd = 287.054; // J/(kg.K)
	double gm = 9.784; // m/s^2
	double g = 9.80665; // m/s^2

	// Meteorological parameters initialization
	double P0, T0, e0, B0, Lam0;
	double Delta_P, Delta_T, Delta_e, Delta_B, Delta_Lam;

	if (Lat <= 15) {
		P0 = 1013.25; T0 = 299.65; e0 = 26.31; B0 = 6.30e-3; Lam0 = 2.77;
		Delta_P = 0; Delta_T = 0; Delta_e = 0; Delta_B = 0; Delta_Lam = 0;
	}
	else if (Lat > 15 && Lat <= 30) {
		P0 = Interpol(Lat, 30, 15, 1017.25, 1013.25);
		T0 = Interpol(Lat, 30, 15, 294.15, 299.65);
		e0 = Interpol(Lat, 30, 15, 21.79, 26.31);
		B0 = Interpol(Lat, 30, 15, 6.05e-3, 6.3e-3);
		Lam0 = Interpol(Lat, 30, 15, 3.15, 2.77);
		Delta_P = Interpol(Lat, 30, 15, -3.75, 0);
		Delta_T = Interpol(Lat, 30, 15, 7, 0);
		Delta_e = Interpol(Lat, 30, 15, 8.85, 0);
		Delta_B = Interpol(Lat, 30, 15, 0.25e-3, 0);
		Delta_Lam = Interpol(Lat, 30, 15, 0.33, 0);
	}
	else if (Lat > 30 && Lat <= 45) {
		P0 = Interpol(Lat, 45, 30, 1015.75, 1017.25);
		T0 = Interpol(Lat, 45, 30, 283.15, 294.15);
		e0 = Interpol(Lat, 45, 30, 11.66, 21.79);
		B0 = Interpol(Lat, 45, 30, 5.58e-3, 6.05e-5);
		Lam0 = Interpol(Lat, 45, 30, 2.57, 3.15);
		Delta_P = Interpol(Lat, 45, 30, -2.25, -3.75);
		Delta_T = Interpol(Lat, 45, 30, 11, 7);
		Delta_e = Interpol(Lat, 45, 30, 7.24, 8.85);
		Delta_B = Interpol(Lat, 45, 30, 0.32e-3, 0.25e-3);
		Delta_Lam = Interpol(Lat, 45, 30, 0.46, 0.33);
	}
	else if (Lat > 45 && Lat <= 60) {
		P0 = Interpol(Lat, 60, 45, 1011.75, 1015.75);
		T0 = Interpol(Lat, 60, 45, 272.15, 283.15);
		e0 = Interpol(Lat, 60, 45, 6.78, 11.66);
		B0 = Interpol(Lat, 60, 45, 5.39e-3, 5.58e-3);
		Lam0 = Interpol(Lat, 60, 45, 1.81, 2.57);
		Delta_P = Interpol(Lat, 60, 45, -1.75, -2.25);
		Delta_T = Interpol(Lat, 60, 45, 15, 11);
		Delta_e = Interpol(Lat, 60, 45, 5.36, 7.24);
		Delta_B = Interpol(Lat, 60, 45, 0.81e-3, 0.32e-3);
		Delta_Lam = Interpol(Lat, 60, 45, 0.74, 0.46);
	}
	else if (Lat > 60 && Lat <= 75) {
		P0 = Interpol(Lat, 60, 75, 1013, 1011.75);
		T0 = Interpol(Lat, 60, 75, 263.65, 272.15);
		e0 = Interpol(Lat, 60, 75, 4.11, 6.78);
		B0 = Interpol(Lat, 60, 75, 4.53e-3, 5.39e-3);
		Lam0 = Interpol(Lat, 60, 75, 1.55, 1.81);
		Delta_P = Interpol(Lat, 60, 75, -0.5, -1.75);
		Delta_T = Interpol(Lat, 60, 75, 14.5, 15);
		Delta_e = Interpol(Lat, 60, 75, 3.39, 5.36);
		Delta_B = Interpol(Lat, 60, 75, 0.62e-3, 0.81e-3);
		Delta_Lam = Interpol(Lat, 60, 75, 0.30, 0.74);
	}
	else if (Lat > 75) {
		P0 = 1013; T0 = 263.65; e0 = 4.11; B0 = 4.53e-3; Lam0 = 1.55;
		Delta_P = -0.5; Delta_T = 14.5; Delta_e = 3.39; Delta_B = 0.62e-3; Delta_Lam = 0.30;
	}

	// Calculate P, T, e, B, Lam using cosine function and DOY
	double P = P0 - Delta_P * cos(2 * PI * (DOY - Dmin) / 365.25);
	double T = T0 - Delta_T * cos(2 * PI * (DOY - Dmin) / 365.25);
	double e = e0 - Delta_e * cos(2 * PI * (DOY - Dmin) / 365.25);
	double B = B0 - Delta_B * cos(2 * PI * (DOY - Dmin) / 365.25);
	double Lam = Lam0 - Delta_Lam * cos(2 * PI * (DOY - Dmin) / 365.25);

	// Calculate Zhyd and Zwet
	double Zhyd = (1e-6 * k1 * Rd * P) / gm;
	double Zwet = (1e-6 * k2 * Rd / (gm * (Lam + 1) - B0 * Rd)) * (e / T);

	// Calculate d_hyd and d_wet
	d_hyd = pow(1 - (B * H / T), g / (Rd * B)) * Zhyd;
	d_wet = pow(1 - (B * H / T), (Lam + 1) * g / (Rd * B) - 1) * Zwet;
}