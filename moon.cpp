#include<iostream>
#include <cmath>
#include <ctime>

using namespace std;

bool isLeapYear(int year);

const int EPOCH_YEAR = 2010;
const double ET_OFFSET_SECONDS = 66.07;

const double YEAR_DURATION_DAYS = 365.242191;

// Sun
const double SOLAR_ECLIPTIC_EPOCH_LONGITUDE = 279.557208;
const double SOLAR_ECLIPTIC_PERIGEE_LONGITUDE = 283.112438;
const double SOLAR_ECCENTRICITY = 0.016705;

// Math
const double PI = 3.141592653589793238463;

// Moon
const double LUNAR_MEAN_EPOCH_LONGITUDE = 91.929336;
const double LUNAR_PERIGEE_EPOCH_LONGITUDE = 130.143076;

double normaliseAngle(double angle)
{
	while (angle < 0 || angle > 360) {
		int direction = (angle < 0) ? 1 : -1;
		angle += direction * 360;
	}
	return angle;
}

int dateToDayOfYear(int day, int month, int year)
{
	if (month <= 2) {
		month -= 1;
		if (isLeapYear(year))
			month *= 62;
		else
			month *= 63;
		month /= 2;
	}	
	else {
		month += 1;
		month = (int) (30.6 * (double) month);
		if (isLeapYear(year))
			month -= 62;
		else
			month -= 63;
	}
	return day + month;

}

int daysSinceEpoch(int day, int month, int year)
{
	int daysToCurrentYear = 0;

	if (year > EPOCH_YEAR) {
		for (int i = EPOCH_YEAR; i != year; i++) {
			daysToCurrentYear += isLeapYear(i) ? 366 : 365;
		}
	}
	if (year < EPOCH_YEAR) {
		for (int i = year; i != EPOCH_YEAR; i++) {
			daysToCurrentYear -= isLeapYear(i) ? 366 : 365;
		}
	}	
	
	return dateToDayOfYear(day, month, year) + daysToCurrentYear;
}

bool isLeapYear(int year)
{
	return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
}

double epochTime(int hour, int minute, int second, int day, int month, int year)
{
	return    ((double) hour / 24.0)
		+ ((double) minute / 1440.0) 
		+ (((double) second + ET_OFFSET_SECONDS ) / 86400.0) 
		+ (double) daysSinceEpoch(day, month, year);
}


double sunN(double epochTime)
{
	double N = (360 / YEAR_DURATION_DAYS) * epochTime;
	return normaliseAngle(N);	
}

double sunMeanAnomaly(double N)
{
	N = N +  SOLAR_ECLIPTIC_EPOCH_LONGITUDE - SOLAR_ECLIPTIC_PERIGEE_LONGITUDE;
	if (N < 0)
		N += 360;
	return N;
}

double sunLongitude(double M, double N) 
{
	double Ec = (360 / PI) * SOLAR_ECCENTRICITY * sin( (M * PI) / 180);
	double l = N + Ec + SOLAR_ECLIPTIC_EPOCH_LONGITUDE;
	
	if (l > 360)
		l -= 360;
	return l;
}

double lunarLongitude(double epochTime)
{
	double l = 13.1763966 * epochTime + LUNAR_MEAN_EPOCH_LONGITUDE;
	return normaliseAngle(l);
}

double lunarAnomaly(double l, double D)
{
	double M = l - 0.111404 * D - LUNAR_PERIGEE_EPOCH_LONGITUDE;
	return normaliseAngle(M);
}

double dtoR(double degrees)
{
	return (PI / 180) * degrees;
}

double evectionCorrection(double l, double L, double m)
{
	double C = l - L;
	return 1.2739 * sin( dtoR(2*C - m) );
}

double annualEquation(double M)
{
	return 0.1858 * sin( dtoR(M) );
}

double thirdCorrection(double M)
{
	return 0.37 * sin( dtoR(M) );
}

double lunarCorrectedAnomaly(double m, double Ev, double Ae, double A3)
{
	return m + Ev - Ae - A3;
}

double centreEquationCorrection(double mm)
{
	return 6.2886 * sin( dtoR(mm) );
}

double fourthCorrection(double mm)
{
	return 0.214 * sin( dtoR(2 * mm) );
}

double lunarCorrectedLongitude(double l, double Ev, double Ec, double Ae, double A4)
{
	return l + Ev + Ec - Ae + A4;
}

double lunarVariation(double ll, double L)
{
	return 0.6583 * sin( dtoR(2 * (ll - L) ) );
}

double lunarOrbitalLongitude(double ll, double V)
{
	return ll + V;
}

double lunarPhase(double lll, double L)
{
	return 0.5 * (1 - cos( dtoR(lll - L) ) );
}


int main()
{
	bool debug = false;
	cout.precision(3);
	
	time_t t = time(0);
	struct tm * now = localtime(&t);
		
	cout << now->tm_mday << "/" << now->tm_mon + 1 << "/" << now->tm_year + 1900 << endl; 

	double D = epochTime(0, 0, 0, now->tm_mday, now->tm_mon + 1, now->tm_year + 1900);
	double N = sunN(D);
	double M = sunMeanAnomaly(N);
	double L = sunLongitude(M, N);
	double l = lunarLongitude(D);
	double m = lunarAnomaly(l, D);
	double Ev = evectionCorrection(l, L, m);
	double Ae = annualEquation(M);
	double A3 = thirdCorrection(M);
	double mm = lunarCorrectedAnomaly(m, Ev, Ae, A3);
	double Ec = centreEquationCorrection(mm);
	double A4 = fourthCorrection(mm);
	double ll = lunarCorrectedLongitude(l, Ev, Ec, Ae, A4);
	double V = lunarVariation(ll, L);
	double lll = lunarOrbitalLongitude(ll, V);
	double F = lunarPhase(lll, L);

	if (debug) {
		cout << "D= " << fixed << D << endl;
		cout << "N= " << fixed << N << endl;
		cout << "M= " << fixed << M << endl;
		cout << "L= " << fixed << L << endl;
		cout << "l= " << fixed << l << endl;
		cout << "m= " << fixed << m << endl;
		cout << "Ev= " << fixed << Ev << endl;
		cout << "Ae= " << fixed << Ae << endl;
		cout << "A3= " << fixed << A3 << endl;
		cout << "mm= " << fixed << mm << endl;
		cout << "Ec= " << fixed << Ec << endl;
		cout << "A4= " << fixed << A4 << endl;
		cout << "ll= " << fixed << ll << endl;
		cout << "V= " << fixed << V << endl;
		cout << "lll= " << fixed << lll << endl;
	}

	cout << "F= " << fixed << F << endl;
	return 0;
}
