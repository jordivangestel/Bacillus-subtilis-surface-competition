// NAME: Jordi van Gestel
// CONTACT: jordivangestel@gmail.com
// DATE: 2022.01.15

// Libraries
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <random>
#include <sstream>
#include <boost/random.hpp>
using namespace std;

// GLOBAL PARAMETERS
int	 S	 = 1;				// random Seed
const int	 T   = 1440;	// total duration of simulation 60*24=1440 minutes
const int	 nx	 = 100;		// total length is 10cmx10cm, so nx=100, would mean that every grid element is 1 mm2.
const double Ri	 = 1000;	// initial resources
const int	 Pi	 = int(0.5*nx/10); // radius of inoculum
const double step_size = 6;	// time step size in seconds at which diffusion is calculated (approximated)

// relevant model parameters; molecule 1 = surfactin, molecule 2 = eps, molecule 3 = blsa
const double p11 = 0;		// privitization of molecule 1 (surfactin) by producer 1
const double p12 = 0;		// privitization of molecule 1 (surfactin) by producer 2
const double p21 = 0.9;		// privitization of molecule 2 (EPS) by producer 1
const double p22 = 0.9;		// privitization of molecule 2 (EPS) by producer 2
const double p31 = 0.1;		// privitization of molecule 3 (BslA) by producer 1
const double p32 = 0.1;		// privitization of molecule 3 (BslA) by producer 2

const double d11 = 0.02;	// production molecule 1 by producer 1
const double d12 = 0.02;	// production molecule 1 by producer 2
const double d21 = 0.002;	// production molecule 2 by producer 1
const double d22 = 0.002;	// production molecule 2 by producer 2
const double d31 = 0.002;	// production molecule 3 by producer 1
const doubled32 = 0.002;	// production molecule 3 by producer 2

const double d 	 = 0.001;	// degradation of molecules (the same for all molecules)
const double c 	 = 1;		// 1 unit of resources consumed per cell doubling
const double m	 = 1.00;	// migration chances of cells in space
const double Th  = 100.0;	// threshold used, before cells can migrate
const double F1	 = pow(10.0,-2.0);		 // Diffusion of molecule 1, dimensionless parameter, diffusion coefficient (mm2/s); reference glucose diffusion 6*10^(-4) mm2/s in 25C water; glycerol in 25C water is 9.3*10^(-4) mm2/s
const double F2	 = pow(10.0,-8.0);		 // diffusion of molecule 2, dimensionless parameter, diffusion coefficient (mm2/s)
const double F3	 = pow(10.0,-4.0);		 // diffusion of molecule 3, dimensionless parameter, diffusion coefficient (mm2/s)
const double FR	 = pow(10.0,-3.0);		 // diffusion of resources, dimensionless parameter, diffusion coefficient (mm2/s)
const double r1	 = exp(log(2.0)/45.0)-1; // growth rate per minute; equates to a cell division every 45 minutes
double		 r2	 = exp(log(2.0)/45.0)-1; // 1.000622*exp(log(2.0)/45.0)-1; // growth rate of surfactin producer

FILE *ParaFile;				// Parameter file
FILE *GridFile;				// Grid information
FILE *SummFile;				// Summary information

// BOOST RANDOM NUMBER GENERATOR /////////////////////////////////////////////////////////////////
typedef boost::mt11213b rndeng;
rndeng rndgen;                      //< The one and only random number generator
int rn(int n)
{
	if(n==1) return 0; //quick escape
	return boost::uniform_int<int>(0,n-1)(rndgen);
}
double ru(){return boost::uniform_01<double>()(rndgen);}
double norm(double sd){return boost::normal_distribution<double>(0,sd)(rndgen);}
void set_seed(unsigned seed){rndgen.seed(rndeng::result_type(seed));}

// NORMAL RANDOM NUMBER GENERATOR FOR BINOMIAL FUNCTION //////////////////////////////////////////
random_device rd;
mt19937 gen(S);
int binomial(int bin_n, double bin_p)
{
	binomial_distribution<> distribution(bin_n,bin_p);
	return distribution(gen);
}

void RV(vector<int> &R, int sizev)
// Randomize vector
{
	vector<int> RN;
	for(int i = 0; i < sizev; ++i) RN.push_back(R[i]);
	random_shuffle(RN.begin(),RN.end());
	R.clear();
	R = RN;
}

struct Grid
{
	// Variables
	vector<vector<double>> M1;		// concentration of molecule 1
	vector<vector<double>> M2;		// concentration of molecule 2
	vector<vector<double>> M3;		// concentration of molecule 3
	vector<vector<int>> P1;			// number of producer 1
	vector<vector<int>> P2;			// number of producer 2
	vector<vector<double>> R;		// resources
	
	void GetN(vector<int> &,vector<int> &, int, int);			// Get neighboring grid elements
	void Init();												// Initialize grid conditions 
};

void Grid::Init()
{
	vector<double>	drow1, drow2;
	vector<int>	irow;
	drow1.assign(nx,0);
	drow2.assign(nx,Ri);
	irow.assign(nx,0);
	M1.assign(nx,drow1);
	M2.assign(nx,drow1);
	M3.assign(nx,drow1);
	P1.assign(nx,irow);
	P2.assign(nx,irow);
	R.assign(nx,drow2);
	
	double r = 1.0 / (2.0 * double(nx) + 1.0);
	double a = sqrt(1.0 + 1.0 / 3.0)*r;
	double h = 0.5*a;
	double xc, yc, xm, ym;
	
	for (int ix = 0; ix < nx; ++ix) for (int iy = 0; iy < nx; ++iy)
	{
		if (iy % 2 == 0) xc = 2.0*r*double(ix) - r;
		if (iy % 2 == 1) xc = 2.0*r*double(ix);
		yc = 1 - a - (a + h)*(double(iy) - 1);
		xm = 2.0*r*double(nx/2.0);
		ym = 1 - a - (a + h)*(double(nx/2.0) - 1);
		
		if (sqrt(pow(xc-xm, 2.0) + pow(yc-ym, 2.0)) < (double(Pi)/double(nx)))
		{
			P1[ix][iy] = 100;
			P2[ix][iy] = 100;
		}
	}
}
Grid G;

//////////////////////////////////////////////////////////////////////////
//							 GRID FUNCTIONS								//
//////////////////////////////////////////////////////////////////////////
void Grid::GetN(vector<int> &gNx, vector<int> &gNy, int gxi, int gyi)
{
	gNx.clear();
	gNy.clear();
	int inits_x[]={gxi-1,gxi+1,gxi,gxi};
	int inits_y[]={gyi,gyi,gyi+1,gyi-1};
	vector<int> gNx_temp(inits_x,inits_x+4);
	vector<int> gNy_temp(inits_y,inits_y+4);
		
	if(gyi%2==0)
	{
		gNx_temp.push_back(gxi+1);
		gNx_temp.push_back(gxi+1);
		gNy_temp.push_back(gyi+1);
		gNy_temp.push_back(gyi-1);
	}
	if(gyi%2==1)
	{
		gNx_temp.push_back(gxi-1);
		gNx_temp.push_back(gxi-1);
		gNy_temp.push_back(gyi+1);
		gNy_temp.push_back(gyi-1);
	}
	assert(gNx_temp.size() == 6 && gNy_temp.size() == 6);
	assert(nx > 1);
	
	vector<int>::iterator Ix = gNx_temp.begin();
	vector<int>::iterator Iy = gNy_temp.begin();
	while(Ix != gNx_temp.end())
	{
		if(*Ix < 0 || *Ix >= nx || *Iy < 0 || *Iy >= nx)
		{

			Ix = gNx_temp.erase(Ix);
			Iy = gNy_temp.erase(Iy);
		}
		else
		{
			++Ix;
			++Iy;
		}
	}
	
	gNx = gNx_temp;
	gNy = gNy_temp;
	gNx_temp.clear();
	gNy_temp.clear();
	
	assert(gNx.size() == gNy.size());
	assert(gNx.size() > 0);
	assert(gNx.size() <= 6);
}

int main()
{
	set_seed(S);
	////////////////////////////
	//		INITIALIZE		  //
	////////////////////////////
	// INITIALIZE GLOBAL (ADMINISTRATION) VECTORS
	G.Init();

	////////////////////////////
	//		RUN PROGRAM		  //
	////////////////////////////
	for(int t = 1; t < (T+1); ++t)
	{
		cout << t << endl;
		vector<vector<double>> tM1, tM2, tM3, tR, dM1, dM2, dM3, dR;
		vector<vector<int>> tP1, tP2, dP1, dP2;
		vector<int> xl, yl, xrl, yrl;
		tM1	 = G.M1;
		tM2	 = G.M2;
		tM3	 = G.M3;
		tP1	 = G.P1;
		tP2	 = G.P2;
		tR	 = G.R;
		
		vector<double>	drow;
		drow.assign(nx,0);
		
		///////////////
		// DIFFUSION //
		///////////////
		for(int s = 0; s < (60/step_size); ++s) // Use Euler's approximate to simulate diffuse for one minute
		{
			dM1.clear(); dM2.clear(); dM3.clear(); dR.clear();
			dM1.assign(nx, drow);
			dM2.assign(nx, drow);
			dM3.assign(nx, drow);
			dR.assign(nx, drow);
			
			#pragma omp parallel for
			for(int xi = 0; xi < nx; ++xi) for(int yi = 0; yi < nx; ++yi)
			{
				vector<int> xl, yl;
				xl.clear(); yl.clear();
				G.GetN(xl,yl,xi,yi);	 // Get all neighboring cells
				for(int j = 0; j < xl.size(); ++j)
				{
					dM1[xi][yi] += 1.0/pow(100.0 / double(nx),2.0) * step_size * F1 * (tM1[xl[j]][yl[j]] - tM1[xi][yi]);
					dM2[xi][yi] += 1.0/pow(100.0 / double(nx),2.0) * step_size * F2 * (tM2[xl[j]][yl[j]] - tM2[xi][yi]);
					dM3[xi][yi] += 1.0/pow(100.0 / double(nx),2.0) * step_size * F3 * (tM3[xl[j]][yl[j]] - tM3[xi][yi]);
					dR[xi][yi]  += 1.0/pow(100.0 / double(nx),2.0) * step_size * FR * ( tR[xl[j]][yl[j]] -  tR[xi][yi]);
				}
			}

			#pragma omp parallel for
			for(int xi = 0; xi < nx; ++xi) for(int yi = 0; yi < nx; ++yi)
			{
				tM1[xi][yi] += dM1[xi][yi];
				tM2[xi][yi] += dM2[xi][yi];
				tM3[xi][yi] += dM3[xi][yi];
				tR[xi][yi]  += dR[xi][yi];
				
				if(tM1[xi][yi] < 0) tM1[xi][yi] = 0;
				if(tM2[xi][yi] < 0) tM2[xi][yi] = 0;
				if(tM3[xi][yi] < 0) tM3[xi][yi] = 0;
				if(tR[xi][yi] < 0)   tR[xi][yi] = 0;
			}
		}
		
		///////////////////////////////////////
		// MOLECULE PRODUCTION & DEGRADATION //
		///////////////////////////////////////
		#pragma omp parallel for
		for(int xi = 0; xi < nx; ++xi) for(int yi = 0; yi < nx; ++yi)
		{
			tM1[xi][yi] = (1 - d)*tM1[xi][yi];
			tM2[xi][yi] = (1 - d)*tM2[xi][yi];
			tM3[xi][yi] = (1 - d)*tM3[xi][yi];
			
			tM1[xi][yi] += d11*double(tP1[xi][yi]) + d12*double(tP2[xi][yi]);
			tM2[xi][yi] += d21*double(tP1[xi][yi]) + d22*double(tP2[xi][yi]);
			tM3[xi][yi] += d31*double(tP1[xi][yi]) + d32*double(tP2[xi][yi]);
		}
		dP1.clear(); dP2.clear();
		dP1 = tP1;
		dP2 = tP2;
		
		///////////////////////////////
		// CELL DIVISION AND SLIDING //
		///////////////////////////////
		for(int xi = 0; xi < nx; ++xi) for(int yi = 0; yi < nx; ++yi)
		{
			int oP1 = binomial(tP1[xi][yi],r1);
			int oP2 = binomial(tP2[xi][yi],r2);
			
			vector<int> V;
			V.clear();
			if(oP1 > 0) for(int pi1 = 0; pi1 < oP1; ++pi1) V.push_back(1);
			if(oP2 > 0) for(int pi2 = 0; pi2 < oP2; ++pi2) V.push_back(2);
			RV(V,V.size());

			for(int iv = 0; iv < V.size(); ++iv)
			{
				if(tR[xi][yi] > c)
				{
					tR[xi][yi]  -= c;
					if(V[iv] == 1)
					{
						if(tM1[xi][yi] > (Th*(1-p11)) && tM2[xi][yi] > (Th*(1-p21)) && tM3[xi][yi] > (Th*(1-p31)) && ru() < m)
						{
							vector<int> xl, yl;
							xl.clear(); yl.clear();
							G.GetN(xl,yl,xi,yi);
							int xr, yr;
							int random_number = rn(xl.size());
							xr = xl[random_number];
							yr = yl[random_number];
							assert(xr <= (xi + 1) & xr >= (xi - 1));
							assert(yr <= (yi + 1) & yr >= (yi - 1));
							dP1[xr][yr] += 1;
						}
						else
						{
							dP1[xi][yi] += 1;
						}
					}
					if(V[iv] == 2)
					{
						if(tM1[xi][yi] > (Th*(1-p12)) && tM2[xi][yi] > (Th*(1-p22)) && tM3[xi][yi] > (Th*(1-p32)) && ru() < m)
						{
							vector<int> xl, yl;
							xl.clear(); yl.clear();
							G.GetN(xl,yl,xi,yi);
							int xr, yr;
							int random_number = rn(xl.size());
							xr = xl[random_number];
							yr = yl[random_number];
							assert(xr <= (xi + 1) & xr >= (xi - 1));
							assert(yr <= (yi + 1) & yr >= (yi - 1));
							dP2[xr][yr] += 1;
						}
						else
						{
							dP2[xi][yi] += 1;
						}
					}
				}
			}
		}
		
		tP1.clear(); tP2.clear();
		tP1 = dP1;
		tP2 = dP2;

		G.M1 = tM1;
		G.M2 = tM2;
		G.M3 = tM3;
		G.P1 = tP1;
		G.P2 = tP2;
		G.R  = tR;
		
	}
	
	return 0;
}


