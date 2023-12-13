
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string>
#include <ctime>
#include <random>  


//Created May 18, 2023
//Updated December 23, 2023

using namespace std;

time_t seconds;

int t, Test, i, j,CurrentDay,LastDay,Extinct;
double Reaction[10000000], CurrentTime, UniSeed,  TimeToRun;
int NumReacts, StartTime;
//Individuals are stored as arrays indexed by age in years
int S[36500],I[36500],M[36500], R[36500], SN[36500], IN[36500],MN[36500], RN[36500], ST,IT,MT,RT,AT,NewCase,WhichDayOfYear;
double TimeStep, T, FS, FI, FM,FR, AST, AIT,AMT, ART;
double b,delta,gamma,omega,lambda,LowLambda,vir,tau,psi,vo;
double Meanb, MeanS, MeanE, MeanI, MeanR;
int ReproDur,MaxAge,SummedOver,BurnYears,LowYears,HighYears,AvgCount,Year;
int Firings[100000000];
double mo, alpha,PropDrop;
double Ibar,IAbar,Mbar,MAbar,CaseRate,Burden;


ofstream out_Densities;
ofstream out_Frequencies;
ofstream out_Means;
ofstream out_Pars;
ofstream out_Burden;
ofstream out_Distri;
ofstream out_TrapDat;
ofstream out_Avgs;

random_device rd;   // non-deterministic generator  
mt19937 gen(rd());  // seed mersenne twister.  

//Define some random number generators
double Big_Rand(double MIN, double MAX)
{
	uniform_real_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Rand_Int(int MIN, int MAX)
{
	uniform_int_distribution<> dist(MIN, MAX); // distribute results between Min and Max
	return dist(gen);
}

int Fish_Rand(double p)
{
	poisson_distribution<> distr(p); // distribute results between 1 and 6 inclusive.  
	return distr(gen);
}





//Define void functions
void Initialize();
void CalculateRates();
void CountReactions();
void ImplementReaction();
void AgeAll();
void SummaryStats();

int main()
{
	//Ask for a filename to which to send simulation output
	string OutputFile;
	cout << "Enter output file name\n";
	cin >> OutputFile;

	//Set the time step for tau leaping
	TimeStep = 1.0;//Currently set at 1 day. Weeks are not fine enough and result in errors. Note that although time is in days, age classes are in years.

	//WARNING: These parameters have units of days unlike the analytical treatment
	//Set static parameters
	b = 500/(365.0);//Set the birth rate
	MaxAge = 100;//Set the maximum possible age in YEARS.
	delta = 1.0/(365.0*60.0);//Set disease independent mortality rate
	gamma = 1.0/14.0;//Set recovery rate
	lambda = 0.4/365.0;//Set the initial force of spillover
	vo = 1.0/365.0;//Set disease dependent mortality rate
	mo = 1.0 / 365.0;//Set the base level rate of transition to the morbid state
	//Ask for the rate of waning immunity
	cout << "Enter omega\n";
	cin>>omega;//Set the rate at which immunity wanes
	LowLambda = 0.05/365.0;//Set the reduced force of spillover after intervention or unintentional anthropogenic change
	//Ask for the paramater tau which defines alpha such that the rate of transition to a diseased state increases by a factor tau between birth and average lifespan
	cout << "Enter tau\n";
	cin>>tau;
	//WARNING: Setting this gets tricky when ages are in years and time units in days. This is done by multiplying alpha by 365 days.
	alpha = 365.0*(delta * mo*(tau - 1));//Set the function defining the rate of transition to the morbid state

	StartTime = 1;//Set the starting time for the simulation

	//Create file names for the various kinds of output all constructed with the root name defined by the user on run
	out_Densities.open("Densities_"+OutputFile);
	out_Burden.open("Burden_" + OutputFile);
	out_Pars.open("Pars_"+OutputFile);
	out_Distri.open("Distribution_" + OutputFile);


	//Set output framework for files
	out_Pars << "b=" << b << ",MaxAge=" << MaxAge << ",delta=" << delta << ",gamma=" << gamma << ",lambda=" << lambda << ",LowLambda=" << LowLambda << ",alpha=" << alpha << "\n";
	out_Pars << "mo=" << mo << ",vo=" << vo <<",omega="<<omega<< "\n";
	out_Burden << "Year,lambda,drop,tau,Ibar,IAbar,Mbar,MAbar,Burden\n";
	out_Densities << "CurrentWeek,b,N,S,I,M,R,NewCases\n";
	out_Distri << "CurrentDay,b,";
	for (i = 0; i < MaxAge; i++)//for the number of ages (this file stores information on the number of individuals in each age and is for debugging purposes only).
	{
		out_Distri <<i<< ",";
	}
	out_Distri << "\n";

	//Flush the output buffers
	out_Densities.flush();
	//Done setting output framework

	Initialize();//Initialize the simulation
	BurnYears = 200;//Set the number of years to run as burn in
	HighYears = 100;//Set the number of years to run with high spillover
	LowYears = 100;//Set the number of years to run with low spillover
	TimeToRun = (BurnYears+LowYears+HighYears) * 365;//Calculate total simulation run time in days
	T = 0;

	//The next couple items are just used to average values over longer periods to generate smoother output
	Burden = 0;
	Mbar = 0;
	Ibar = 0;
	AvgCount = 0;
	Year = 0;
	do//start time loop for simulation
	{
		//Averager (start this after burn in)
		if (T >= BurnYears*365)
		{
			CaseRate = NewCase / (1.0*AT);//Divide the total number of disease cases in this time step (NewCase) by the total population size
			Burden = Burden + CaseRate;//Sum these per capita rates up within the year to get annual burden
			Mbar = Mbar + FM;
			Ibar = Ibar + FI;
			MAbar = MAbar + AMT;
			IAbar = IAbar + AIT;
			AvgCount++;

			if (AvgCount == 365)//once a year has passed, average things up and reset for the next year
			{
				Mbar = Mbar / 365.0;
				MAbar = MAbar / 365.0;
				Ibar = Ibar / 365.0;
				IAbar = IAbar / 365.0;
				out_Burden << Year << "," << lambda << "," << PropDrop << "," << tau << "," << 100 * Ibar << "," << IAbar << "," << 100 * Mbar << "," << MAbar << "," << Burden << "\n";//Note that output units on age here have been changed to years...
				AvgCount = 0;
				Mbar = 0;
				MAbar = 0;
				Ibar = 0;
				IAbar = 0;
				Burden = 0;//Clear the burden counter out so it is ready for the next year
				Year++;
			}
		}
		//End of annual averager
		
		CurrentDay = int(T);//Discretize time
		WhichDayOfYear = CurrentDay % 365;//Record day of the year

		cout << TimeToRun - T << "\n";
		//Execute Reactions
		CalculateRates();
		CountReactions();
		ImplementReaction();
		SummaryStats();
		

		//Implement spillover reduction at the defined time by changing lambda
		if (T > (HighYears+BurnYears) * 365)
		{
			lambda = LowLambda;
		}


		//Output current state of the population before aging them
		out_Densities << T << "," << b << "," << AT  << "," << ST  << "," << IT << "," << MT << "," << RT <<","<<NewCase<<"\n";

		if (CurrentDay % 365 == 0)//output full age distribution once a year and age everyone once per year
		{
			out_Distri << T << "," << b << ",";
			for (i = 0; i < MaxAge; i++)//for all ages
			{
				out_Distri << S[i] << ",";
			}
			out_Distri << "\n";

			AgeAll();//age everyone
		}

		T = T + TimeStep; //Increment the time step


	} while (T < TimeToRun);//End time loop when the predfined end-time is reached

}

void Initialize()//This function initializes the population
{

	//Clear the arrays
	for (i = 0; i < MaxAge; i++)
	{
		S[i] = 0;
		I[i] = 0;
		M[i] = 0;
		R[i] = 0;
	}
	//Start the simulation with b newborn, susceptible individuals
	S[0] = b;
	I[0] = 0;
	M[0] = 0;
	R[0] = 0;

}

void CalculateRates()//This function calculates the rate of each reaction type
{
	int CR;

	CR = 0;//CR is just a reaction counter

	//births occur
	Reaction[CR] = b ;//A birth occurs
	CR++;

	//Disease independent death
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = delta * S[i];
		CR++;
		Reaction[CR] = delta * I[i];
		CR++;
		Reaction[CR] = delta * M[i];
		CR++;
		Reaction[CR] = delta * R[i];
		CR++;
	}

	//Disease dependent death
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = vo * M[i] ;
		CR++;
	}

	//Spillover
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = lambda * S[i];
		CR++;
	}

	//Pathogen recovery for infected class
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = gamma * I[i];
		CR++;
	}

	//Pathogen recovery for morbid class
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = gamma * M[i];
		CR++;
	}

	//Progression to disease
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = (mo+alpha*i) * I[i];
		CR++;
	}

	//Waning immunity
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		Reaction[CR] = omega* R[i];
		CR++;
	}

	NumReacts = CR;

	if (NumReacts > 10000000)//Check to make sure we have not overflowed the reaction rate array
	{
		cout << "WARNING: Reaction Overflow\n";
		cin >> Test;
	}

	//Error checks follow
	for (i = 0; i < MaxAge; i++)
	{
		if (S[i] < 0 || I[i] < 0 || M[i]<0||R[i] < 0)//If any density goes negative, throw an error and await acknowledgement
		{
			cout << "WARNING: Negative densities have occured\n";
			cin >> Test;
		}
	}
}

void CountReactions()//Draw the number of reactions of each type that occur at random
{
	double p;
	int CR;

	//Step through all possible reactions and draw the number that occur at random from a poisson with mean rate*t
	for (CR = 0; CR < NumReacts; CR++)
	{
		p = Reaction[CR] * TimeStep;
		if (p > 0)
		{
			Firings[CR] = Fish_Rand(p);
		}
		else
		{
			Firings[CR] = 0;
		}

	}
}


void ImplementReaction()//Update the abundance of each type
{
	int CR;

	CR = 0;//CR is just a reaction counter
	NewCase = 0;//Reset the counter for new cases this week
	
	while (Firings[CR] > 0)//Birth
	{
		S[0]++;
		Firings[CR]--;
	}
	CR++;

	for (i = 0; i < MaxAge; i++)//Disease independent death 	
	{
		while (Firings[CR] > 0 && S[i] > 0)
		{
			S[i]--;
			Firings[CR]--;
		}
		CR++;
		while (Firings[CR] > 0 && I[i] > 0)
		{
			I[i]--;
			Firings[CR]--;
		}
		CR++;
		while (Firings[CR] > 0 && M[i] > 0)
		{
			M[i]--;
			Firings[CR]--;
		}
		CR++;
		while (Firings[CR] > 0 && R[i] > 0)
		{
			R[i]--;
			Firings[CR]--;
		}
		CR++;
	}

	for (i = 0; i < MaxAge; i++)//Disease dependent death
	{
		while (Firings[CR] > 0 && M[i] > 0)
		{
			M[i]--;
			Firings[CR]--;
		}
		CR++;
	}
	
	for (i = 0; i < MaxAge; i++)//Spillover
	{
		while (Firings[CR] > 0 && S[i] > 0)
		{
			S[i]--;
			I[i]++;
			Firings[CR]--;
		}
		CR++;
	
	}
	for (i = 0; i < MaxAge; i++)//Pathogen recovery for the I class
	{
		while (Firings[CR] > 0 && I[i] > 0)
		{
			I[i]--;
			R[i]++;
			Firings[CR]--;
		}
		CR++;
	}
	for (i = 0; i < MaxAge; i++)//Pathogen recovery for the M class
	{
		while (Firings[CR] > 0 && M[i] > 0)
		{
			M[i]--;
			R[i]++;
			Firings[CR]--;
		}
		CR++;
	}
	//Progression to disease 
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		while (Firings[CR] > 0 && I[i] > 0)
		{
			NewCase++;
			I[i]--;
			M[i]++;
			Firings[CR]--;
		}
		CR++;
	}
	//Waning immunity 
	for (i = 0; i < MaxAge; i++)//for the number of age classes
	{
		while (Firings[CR] > 0 && R[i] > 0)
		{
			R[i]--;
			S[i]++;
			Firings[CR]--;
		}
		CR++;
	}
}

void AgeAll() // This is currently set such that all individuals of MaxAge are lost from the system implying certain death at MaxAge. Note that this can lead to discrepancies with analytical solutions due to their assuming an exponential distribution of age.
{
	for (i = 1; i < MaxAge; i++)//Advance everyone into temporary arrays
	{
		SN[i] = S[i - 1];
		IN[i] = I[i - 1];
		MN[i] = M[i - 1];
		RN[i] = R[i - 1];
	}
	S[0] = 0;
	I[0] = 0;
	M[0] = 0;
	R[0] = 0;
	for (i = 1; i < MaxAge; i++)//Reset arrays
	{
		S[i] = SN[i];
		I[i] = IN[i];
		M[i] = MN[i];
		R[i] = RN[i];
	}
}

void SummaryStats()
{
	//First, summarize densities/abundances of each class
	ST = 0;
	IT = 0;
	MT = 0;
	RT = 0;
	for (i = 0; i < MaxAge; i++)
	{
		ST = ST + S[i];
		IT = IT + I[i];
		MT = MT + M[i];
		RT = RT + R[i];
	}
	AT = ST+IT+MT+RT;


	//Next, summarize the frequencies of each class
	FS = ST / (1.0*AT);
	FI = IT / (1.0*AT);
	FM = MT / (1.0*AT);
	FR = RT / (1.0*AT);

	//Summarize the average age of each class
	AST = 0;
	AIT = 0;
	AMT = 0;
	ART = 0;
	for (i = 0; i < MaxAge; i++)
	{
		if (ST > 0)
		{
			AST = AST + i * (S[i] / (1.0*ST));
		}
		else
		{
			AST = -1;
		}
		if (IT > 0)
		{
			AIT = AIT + i * (I[i] / (1.0*IT));
		}
		else
		{
			AIT = -1;
		}
		if (MT > 0)
		{
			AMT = AMT + i * (M[i] / (1.0*MT));
		}
		else
		{
			AMT = -1;
		}
		if (RT > 0)
		{
			ART = ART + i * (R[i] / (1.0*RT));
		}
		else
		{
			ART = -1;
		}
	}
	
}









