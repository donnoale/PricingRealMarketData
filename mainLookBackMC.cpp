// By Alessandro Donno, 2017

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <string>
#include <chrono>
#include <ctime>
#include <windows.h>
#include <mmsystem.h>

using namespace std;

double riskyStockPoint(double S0, double sigma, double T, double r)
{
	
// this function computes the price of a stock after a given time period in a standard GBM framework

unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();

mt19937 e(seed);	// A Mersenne Twister pseudo-random generator

normal_distribution<double> distrN(0.0, 1.0);		// Random number distribution that produces floating-point values according to a normal distribution
double number = distrN(e);
double result = S0 * exp((r - (sigma*sigma*0.5)) * T + sigma * pow(T, 0.5) * number);
return result;
}


int main() {

// this function runs the Monte Carlo simulation of a lookback option following a standard GBM process

string optionType;
double K;
double S0;
double sigma;
double T;
double r;
const int b = 501;
const int c = 100001;
double S[b];
double payoff[c];
double discPayoff[c];
double dt;
int n;
int nIter;
double price;
double radius;
double lcl;
double ucl;
double lcl2;
double sampleStdDev;
double sumdiff;

cout<<"Type option flavor (C/P)"<<endl;
cin>>optionType;
cout<<"Type S0"<<endl;
cin>>S0;
cout<<"Type sigma"<<endl;
cin>>sigma;
cout<<"Type T"<<endl;
cin>>T;
cout<<"Type r"<<endl;
cin>>r;
cout<<"Type number of time steps"<<endl;
cin>>n;
cout<<"Type number of iterations"<<endl;
cin>>nIter;
cout<<endl;

// start the timer
clock_t begin = clock();

if (optionType=="C" || optionType=="c" || optionType=="P" || optionType=="p")	{

S[0]=S0;
dt=T/static_cast<double>(n);
discPayoff[0]=0.0;

for (int a=1; a<=nIter; a++)
{
	for (int i=1; i<=n; i++)
	{
		S[i] = riskyStockPoint(S[i-1], sigma, dt, r);
	}
	if (optionType=="C" || optionType=="c")
		{
			K = *min_element(S,S+n+1);		// Returns an iterator pointing to the element with the samllest value in a vector
			payoff[a] = max(0.0, S[n] - K);
		}
		else
		{
			K = *max_element(S,S+n+1);		// Returns an iterator pointing to the element with the largest value in the vector
			payoff[a] = max(0.0, K - S[n]);
		}
	
	discPayoff[a] = discPayoff[a-1] + (exp(-r * T) * payoff[a]);

}

price = discPayoff[nIter] / static_cast<double>(nIter);
}
else	{
			cout<<"Unknown option type: "<<optionType<<endl<<"Type C to get the value of a Call Option"<<endl<<"Type P to get the value of a Put Option"<<endl<<endl;
            price = 0.0;
}
	
cout<<"The estimated price of this option is:"<<endl;
cout<<price<<endl<<endl;

for (int a=1; a<=nIter; a++)
{
	sumdiff = sumdiff + pow(((exp(-r * T) * payoff[a]) - price), 2.0);
}

sampleStdDev = pow((1.0 / (static_cast<double>(nIter) - 1.0)) * sumdiff, 0.5);
radius = (sampleStdDev / pow(static_cast<double>(nIter), 0.5)) * 1.96;
cout<<"The radius of the 95% confidence interval is:"<<endl;
cout<<radius<<endl<<endl;

lcl = price - radius;
ucl = price + radius;
cout<<"The 95% confidence interval is:"<<endl;
cout<<lcl << " < price < " << ucl<<endl<<endl;

clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
cout<<"This code ran successfully in "<<elapsed_secs<<" seconds"<<endl<<endl;

system("pause");

return 0;
}
