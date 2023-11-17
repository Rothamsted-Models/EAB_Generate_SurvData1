#include "EABGrid.h"
#include "EABModel.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <direct.h>
#include <list>

#include "D:\\Users\\milne\\AppData\\Local\\NAG\\NL29\\nlw3229del\\include\\nag.h"


//Pop dynamics EAB params
// Time
double P_LarvalDays = 56; // 8 weeks of development according to Duan et al 2013

// Larvae
double P_LarvToAdultStart = 0.25; // Grid::GetPropOverWinterLarvaeStart();//LarvaeToAdultInOneYear(Lyons, 2015)
double P_Cannibalism = 0.0001986334; //  0.00025; // 0.0075 Death rate due to cannibalism per m ^ 2 of phloem surface area(commented ones are obtained from the transformations in PPT)
// (Duan et al 2013 fig2)  sin ^ 2(0.000578) = 0.0000003341 Green ash; sin ^ 2(0.000404) = 0.000000163216 Tropical ash
double P_Starvation =  0.0036; // 0.00035; // 0.0095;% Death rate due to starvation 0.0036 is US and represents worst case scenario. UK may be bigger could be  0.8/56 Shouwalter(?)
// (Duan et al 2013 fig1 and fig2) sin ^ 2(0.000220) = 0.00000004834 Green ash, sin ^ 2(0.000105) = 0.00000001102 Tropical ash
double LarvalDeath1 = exp(-P_Starvation * P_LarvalDays);
double LarvalDeath2 = (P_Cannibalism / P_Starvation) * (1 - exp(-P_Starvation * P_LarvalDays));




// Adults
double P_AdultDeath = 0.17 / 0.279; //ProbDeathB4Reproduction 0.17
double P_NumberEggs = 26.58; //- 35 LarvaeOffspring - female depostis eggs in bark cracks in small groups of two or 3 eggs(Mcquarrie et al 2015)
double P_LayingEgg = 0.0001; // This is the probability of laying an egg range [0.0005 - 0.00004]

double treekappa = 0.0138939;// this is now per emergence hole per m2 //0.25; // 0.25; //tree death rate once infected Put this to 1 to avoid tree death

/////////////////////////////////////////////////////////////////////////////////////////////////

void InitialiseModel( )
{
	
	Grid::ClearData();
	FlightP::ClearData();
	//Initialises some outputfiles that we shall write our results in. 
	//Currently commented out as we don't know what we want to write out
	
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	strcat_s(Npath_buffer, "\\OutFiles");
	int  dirMade=_mkdir(Npath_buffer);
	
	SetLand(0); //0 = read from files		
	
	
	

	//In SetLand - we have initialised some trees in a cell as one year infected so we need to pop the EAB on that. This is what the 
	//code below is doing - don't panic its fine
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			if ((Grid::GetCrop(icount, jcount,4)+ Grid::GetCrop(icount, jcount, 7)) >0) ////Status of trees. Status 4 is 'removed which includes never planted there. If this is the case there will be no bug
			{

				
				Grid::SetNumG(icount, jcount, 10.0); //Sets number healthy bugs
				Grid::SetOWLav(icount, jcount, 0);
				Grid::SetTempBug(icount, jcount, 0);
				Grid::SetPLarv(icount, jcount, P_LarvToAdultStart);
				Grid::SetInvasion(icount, jcount, -1);
			
				
			}
			else
			{
				Grid::SetNumG(icount, jcount,  0); //Sets number healthy bugs
				Grid::SetOWLav(icount, jcount, 0);
				Grid::SetTempBug(icount, jcount, 0);
				Grid::SetPLarv(icount, jcount, P_LarvToAdultStart);
				Grid::SetInvasion(icount, jcount, -1);
				
			}
			
		
		}
	}

	
	FlightP::Initialise(); // This initialises the expotential decay for flight
	
	



}

void DeleteModel()
{
	Grid::ClearData();
	FlightP::ClearData();
	

}

void SetLand(int readin)
{
	int iflag(0);
	Grid::SetLand(readin, iflag);
	if (iflag==1)
		throw std::logic_error("Error in Planta assignment");
	else if (iflag==2)
		throw std::logic_error("Error creating Planta landscape");
	else if (iflag==3)
		throw std::logic_error("Error creating Block assignment");
	
}

//Eggs are laid on tips of growing shoots on and between unfurling leaves. Females lay 300 to 800 eggs during their lifetime. 
//Nymphs pass through five instars. The total life cycle requires from 15 to 47 days, depending on environmental factors such as temperature and season. 
//The adults may live for more than a month. There is no diapause, but populations are typically low in the winter or during dry periods. 
//There are 9 to 10 generations a year, with up to 16 observed under observation in field cages. 


//Note this is an annual time step not monthly as the name suggests!
void ModelMonth(int imth, double randInf, double kill_rate)
{

	//This runs the bugs life cycle
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			Generation(icount, jcount, imth, kill_rate) ; //grows bugs and calculates the number of ingected ones based on tree health
		}
	}

//	WriteBugs(imth); // this function writes a map of infected bugs for total see Tbugs

	Flight(); //Bugs distribute
	

	for (int icount = 0; icount < Grid::GetNumRows(); icount++)
	{
		for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
		{
			int FellingType(0);
			double propFelled(0);
			Grid::GetTreeFellingStats(icount, jcount, FellingType, propFelled);
			if (FellingType != 0)
				int junk = 0;
			AdjustEABStatusPostFlight(icount, jcount, FellingType, propFelled);
		}
	}




	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
			
			Tree_Generation(icount, jcount, imth); //ADB Disease - spread not modelled yet
			
		}
	}
				
}

void Generation(int irow, int jcol, int imth, double kill_rate) // Model the population dynamics over a month 
{
	
	
	double TAdult = 0;
	double TLarvae = 0;
	double P_LarvToAdult(0);

	//Get the adults from last year 
	TAdult = Grid::GetNumG(irow, jcol); // Number of adults in cell
	//Get the larvae from last year
	Grid::GetOWLav(irow, jcol, TLarvae); // Density larvae per tree NOTE DIFFERENCE IN UNITS BETWEEN ADULTS AND LARVAE!
	

	if ((TAdult + TLarvae) > 0)
	{
		Grid::GetPLarv(irow, jcol, P_LarvToAdult);
		//Increase in bugs per time step
		double spray = 0.0; // If they buy into control strategy and its the month when control is applied then spray
		
		double tempTreeStates[12];
		for (int icount = 0; icount < Grid::GetMaxNumStatus(); icount++)
		{
			tempTreeStates[icount] = Grid::GetCrop(irow, jcol, icount); // These are proportions
		}

		//double k_cap = Grid::GetTreeDensity(irow, jcol);
		double Insecticide = spray * kill_rate; // control applied to current time step, maybe needeed later to kill bugs
		
		
		//Get Infected trees 
		double propENoEAB = tempTreeStates[0] + tempTreeStates[1] + tempTreeStates[2] + tempTreeStates[3];
		if (propENoEAB > 1)
			propENoEAB = 1;
		double propEInfGTYr = tempTreeStates[8] + tempTreeStates[9] + tempTreeStates[10] + tempTreeStates[11]; // GT a year of EAB
		double NumTreesHealthy = Grid::GetTreeDensity(irow, jcol)* propENoEAB;  // Note GetTreeDensity returns tree number!
		double NumTreesInfested = Grid::GetTreeDensity(irow, jcol) * (1-propENoEAB);
		double NumTreesInfGTYr = Grid::GetTreeDensity(irow, jcol) * propEInfGTYr;
		double phloemPerTree=Grid::GetPhloemPerTree();
		double AdultsM2Phloem = TAdult / (NumTreesInfested* phloemPerTree); //The is the average across all trees newly infected and old
		
		//Adults that emerge from trees infected for a while
		double AdultsM2Phloem_2 = (((1 - P_LarvToAdult) * (1 - P_AdultDeath) * P_NumberEggs * LarvalDeath1) * AdultsM2Phloem + LarvalDeath1 * TLarvae) / (1 + (LarvalDeath2 * (1 - P_AdultDeath) * P_NumberEggs * AdultsM2Phloem) + (LarvalDeath2 * TLarvae));
		//Adults that emerge from trees infected only last year
		double AdultsM2Phloem_1 = (((1 - P_LarvToAdult) * (1 - P_AdultDeath) * P_NumberEggs * LarvalDeath1) * AdultsM2Phloem) / (1 + (LarvalDeath2 * (1 - P_AdultDeath) * P_NumberEggs * AdultsM2Phloem));
		
		//This is Larvae per m2 phloem 
		TLarvae = ((P_LarvToAdult) * (1 - P_AdultDeath) * P_NumberEggs * LarvalDeath1 * AdultsM2Phloem) / (1 + (LarvalDeath2 * (1 - P_AdultDeath) * P_NumberEggs * AdultsM2Phloem) + LarvalDeath2 * TLarvae);

		//Then the total number of adults emerging from this cell is 
		//Delta change in trees(NumTreesInfested - NumTreesInfsLT 
		double Delta = NumTreesInfested - NumTreesInfGTYr; // Total number infested - number infested for a long time

		if (Delta < 0)
		{
			if (Delta > -0.00000000001) // rounding issue ... need to dig deeper here as 
				Delta = 0;
			else
				throw std::logic_error("Error in ModelMonth");

		}
			

		TAdult = (Delta * AdultsM2Phloem_1 + NumTreesInfGTYr * AdultsM2Phloem_2)* phloemPerTree;

		P_LarvToAdult = -0.194 * log(TLarvae/ P_LarvToAdult) + 0.9862;

		if ((P_LarvToAdult <= 0.05))
		{
			P_LarvToAdult = 0.05;
		}
		else if ((P_LarvToAdult >0.95))
		{
			P_LarvToAdult = 0.95;
		}

		Grid::SetPLarv(irow, jcol, P_LarvToAdult); //ALICE PUT BACK
		Grid::SetNumG(irow, jcol, TAdult);
		Grid::SetOWLav(irow, jcol, TLarvae);
		Grid::SetEmergenceDensity(irow, jcol, TAdult / (NumTreesInfested * phloemPerTree));

			
			
	}
	else
	{
		Grid::SetNumG(irow, jcol, TAdult);
		Grid::SetOWLav(irow, jcol, TLarvae);		
	}

	//ofm.close();
}

void Tree_Generation(int irow, int jcol, int iyr) //progesses tree health status in a time step
{


	//This is for ADB

	//char path_buffer[_MAX_PATH];
	//char Npath_buffer[_MAX_PATH];
	//_getcwd( path_buffer, _MAX_PATH );
	//strcpy_s(Npath_buffer, path_buffer);
	//strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");

	// three states of EAB - none, present not detected, detected	
	double H_trees[3] = { Grid::GetCrop(irow, jcol, 0), Grid::GetCrop(irow, jcol, 4), Grid::GetCrop(irow, jcol, 8) }; //Healthy trees
	double E_trees[3] = { Grid::GetCrop(irow, jcol, 1), Grid::GetCrop(irow, jcol, 5), Grid::GetCrop(irow, jcol, 9) }; // Infected not infectious - no EAB
	double I_trees[3] = {Grid::GetCrop(irow, jcol, 2), Grid::GetCrop(irow, jcol, 6), Grid::GetCrop(irow, jcol, 10) }; // Infectious not symptomatic - no EAB
	double S_trees[3] = { Grid::GetCrop(irow, jcol, 3), Grid::GetCrop(irow, jcol, 7), Grid::GetCrop(irow, jcol, 11) }; // Symptomatic - no EAB


	//This needs to be some infection spread eventually
	double Bug_I = 0.01; 

	
}

void AdjustEABStatusPostFlight(int irow, int jcol, int FellingType, double propFelled) // Felling type:  0 = no felling; 1 = fell EAB infected; 2= sanitation felling; propFelled tells the proportion of ill trees to fell
{
	//###########################################################################
			//Re-adjust tree infection status. We will want to think about this - for now I'm moving last years EAB undetected into detected.
			//And putting any new landings into present undetected. This will need to change

			//Consider those trees with ash die back: // This is either Healthy=0, 
			//No EAB +Exposed  by ADB but not infectious = 1, 
			//No EAB +infected by ADB infectctious no symtoms=2, 
			//No EAB + Symptoms by ADB =3, 
				//EAB but only in last year + no ADB =4, 
				//EAB but only in last year +Exposed by ADB but not infectious = 5, 
				//EAB but only in last year +infected by ADB infectctious no symtoms=6, 
				// EAB but only in last year + Symptoms by ADB =7, 
					//EAB more than a year + no ADB =8, 
					//EAB more than a year +Exposed by ADB but not infectious = 9, 
					//EAB  more than a year +infected by ADB infectctious no symtoms=10, 
					// EABmore than a year + Symptoms by ADB =11, 
	double numTrees = Grid::GetTreeDensity(irow, jcol);
	
	
	if (numTrees > 0)
	{
		if (FellingType == 2) // kill all trees
		{
			double TotalSurvivingTrees = 0;
			Grid::SetCrop(irow, jcol, 0, 1); // set proportions to all healthy - we have no trees left! 
			for (int iicount = 1; iicount < Grid::GetMaxNumStatus(); iicount++)
			{
				Grid::SetCrop(irow, jcol, iicount, 0);
			}
			Grid::SetTreeDensity(irow, jcol, TotalSurvivingTrees);

		}
		else
		{

			double TAdult = Grid::GetNumG(irow, jcol); // this is after flight not emerge
			double TreeDeath = Grid::GetEmergenceDensity(irow, jcol) * treekappa; // NB for per emergence hole death rate
			//Get the eggs laid by healthy females 
			if (TAdult > 0)
			{

				double tempTreeStates[12];
				for (int icount = 0; icount < Grid::GetMaxNumStatus(); icount++)
				{
					tempTreeStates[icount] = Grid::GetCrop(irow, jcol, icount);
				}

				//Add up last times EAB status 
				//This is tricky as proportions as we will loose some trees now due to death 
				//Lets run over the four states of ADB 
				double EnoEAB[4] = { 0,0,0,0 };
				double NewInfEAB[4] = { 0,0,0,0 };
				double OldInfEAB[4] = { 0,0,0,0 };
				double Eaten[4] = { 0,0,0,0 };
				double TotEABtrees = 0.0; // new NB

				for (int iADB = 0; iADB < 4; iADB++)
				{
					// Uninfected to infected state 
					EnoEAB[iADB] = tempTreeStates[iADB] * numTrees * exp(-P_LayingEgg * TAdult); //trees left uninfested
					NewInfEAB[iADB] = tempTreeStates[iADB] * numTrees - EnoEAB[iADB]; // trees just infested = difference between then and now 
					//test to make sure that is positive 
					if (NewInfEAB[iADB] < 0)
					{
						NewInfEAB[iADB] = 0;
					}
					OldInfEAB[iADB] = (tempTreeStates[iADB + 4] + tempTreeStates[iADB + 8]) * numTrees; // number of trees infested for over a year
					if (FellingType == 1) // we fell a propotion of sick trees
					{
						OldInfEAB[iADB] = OldInfEAB[iADB] * propFelled;
					}
					//Eaten[iADB] = OldInfEAB[iADB] * (1 - treekappa);  // NB comment out
					if (TreeDeath > 0) {
						Eaten[iADB] = OldInfEAB[iADB] * TreeDeath;  // NB swap to this
						OldInfEAB[iADB] = OldInfEAB[iADB] - Eaten[iADB]; // adjust old infected trees for those that have dies // NB comment out
					}

				}

				double TotalSurvivingTrees = 0.0;
				for (int kcount = 0; kcount < 4; kcount++)
				{


					TotalSurvivingTrees = TotalSurvivingTrees + EnoEAB[kcount] + NewInfEAB[kcount] + OldInfEAB[kcount];
				}

				if (TotalSurvivingTrees > numTrees + 0.0000000001)
				{
					throw std::logic_error("How did that happen ");
				}


				if (TotalSurvivingTrees <= 0)
				{
					TotalSurvivingTrees = 0;
					Grid::SetCrop(irow, jcol, 0, 1); // set proportions to all healthy - we have no trees left! 
					for (int iicount = 1; iicount < Grid::GetMaxNumStatus(); iicount++)
					{
						Grid::SetCrop(irow, jcol, iicount, 0);
					}


				}
				else
				{

					//First the proportion of trees that had EAB for just a year now gets added to the trees that had it before. This is last years I_n becomine I_n-1
					for (int kcount = 0; kcount < 4; kcount++)
					{
						Grid::SetCrop(irow, jcol, 0 + kcount, EnoEAB[kcount] / TotalSurvivingTrees);
						Grid::SetCrop(irow, jcol, 4 + kcount, NewInfEAB[kcount] / TotalSurvivingTrees);
						Grid::SetCrop(irow, jcol, 8 + kcount, OldInfEAB[kcount] / TotalSurvivingTrees);
					}


				}

				Grid::SetTreeDensity(irow, jcol, TotalSurvivingTrees);
			}

		}

		



	}
}

void Flight()
{
	//Disperse first
	//Clear Temp Vector
	Grid::ClearTempBugs();
	int season(1); //summer
	int post(0); //pre=0 post=1
	
	
	BugsJump(); // first see if we have any jumpers
		
//	WriteBugsT(1);
	
	//Bugs Fly 
	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
				
				double I_Bug=Grid::GetNumG(icount, jcount); 
				if (I_Bug>0)
					Distribute2(icount, jcount, I_Bug);

				//Do we get a jump 
				//1.	Probability one jumps on a car = Function of number of EAB
				//	2.	Distance travelled – sampled from Nathan’s graphic
				//	3.	Direction of travel(assume all equally likely)

		}
	}
	//Set grid of Bugs with new values 


	for (int icount=0; icount<Grid::GetNumRows(); icount++)
	{
		for (int jcount=0; jcount<Grid::GetNumCols(); jcount++)
		{
				Grid::SetNumG(icount, jcount,  Grid::GetTempBug(icount, jcount));  //infected
			//	Grid::SetOWLav(icount, jcount, Grid::GetNumG(icount, jcount));
		}
	}
	
	
	//WriteBugsT(2);
	
}

void BugsJump()
{
	Grid::ClearTempBugs();
	//Probability 
// there is something i dont get about how this works!!
	// Question for Alice - should we swap to Mersene Twister genid =3 (this is R standard now), 
	// and do we need to make sure we only intitalise once, otherwise sequence not neccesarily random
	// i have removed unnessary steps here, but there are multiple initialisations in different functions
	// i have added extra nag headers to EABGrid.h - to cover rand poisson, 
	double urand[3]; //first is probability that we set a bug, second and third define which cell in which plant we choose
	
	int ifail = 1;
	int lstate(0);
	int subid(1);
	int* state;
	state = new int[lstate];
	int Lseed = 1;
	int seed[1];
	seed[0] = 1;
	G05KGF(3, subid, state, lstate, ifail);
	//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
	delete[] state;
	state = new int[lstate];
	ifail = 1;
	G05KGF(3, subid, state, lstate, ifail);
	//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);



	// Poisson
	int nJumps[1];
	double *r = 0;
	int lr; 
	NagError fail;
	
	//Bugs Fly 
	for (int icount = 0; icount < Grid::GetNumRows(); icount++)
	{
		for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
		{
			
			double I_Bug = Grid::GetNumG(icount, jcount); // Get number of adults 
			if (I_Bug > 1) // there must be at least one bug
			{
				double Pmean = 0.0002 * I_Bug; // needs lamda to set up generator memory
				lr = 30 + 20 * sqrt(Pmean) + Pmean;
				r = NAG_ALLOC(lr, double);
				INIT_FAIL(fail);
				g05tjc(Nag_InitializeAndGenerate, 1, Pmean, r, lr, state, nJumps, &fail); // poisson number of jumps
				while (nJumps[0] > 0) { // if jumps happen
					G05SAF(3, state, urand, ifail);
					if (ifail > 0)
						throw std::logic_error("Error Opps");
						//Sample Distance
					//options
						//Calc inverse of CDF of exponential (because its easy) 
						double myDistance = -  47 * log(1 - urand[1]); // Cambridge 20km year for 1997 to 2002; then 47km yeardouble 
						//double myDistance = -  87 * log(1 - urand[1]); // nb changes
						//double myDistance = -log(urand[1])/0.06; // USA
						//double myDistance = 10.125 * ((2 * urand[0] - 1) / urand[1]);//cauchy Russia
						//Sample theta 
						double theta = 2 * 3.141572 * urand[2];
						double myGridLen = Grid::GetGridLen();
						int myX = floor(myDistance * sin(theta) / myGridLen) + jcount;
						int myY = floor(myDistance * cos(theta) / myGridLen) + icount;

						if (myX >= Grid::GetNumCols())
							myX = Grid::GetNumCols() - 1;
						if (myY >= Grid::GetNumRows())
							myY = Grid::GetNumRows() - 1;
						if (myX < 0)
							myX = 0;
						if (myY < 0)
							myY = 0;

						//Look to see if there is a bug and add one 
						double myCrop = Grid::GetTreeDensity(myY, myX); //get treedensity
						int numMove = 50;
						if (I_Bug < numMove)
							numMove = 1;// if not enough bugs then just jump one, maybe add error?

						if (myCrop > 0.0 ) //Trees to land on 
						{
							Grid::SetNumG(icount, jcount, I_Bug - numMove); // add catch for incase not 50 beetles there
							double J_Bug = Grid::GetNumG(myY, myX) + numMove;
							//Grid::SetNumG(myY, myX, J_Bug);
							Grid::SetTempBug(myY, myX, J_Bug);
						}
						else
						{ // edit to stop EAB travelling to the beach
							nJumps[0] = nJumps[0] + 1; // if no crop start again - cant jump to sea
						}
					nJumps[0] = nJumps[0] - 1;
				}
				NAG_FREE(r);


			}
			
			//Do we get a jump 
			//1.	Probability one jumps on a car = Function of number of EAB
			//	2.	Distance travelled – sampled from Nathan’s graphic
			//	3.	Direction of travel(assume all equally likely)

		}
	}
	//Set grid of Bugs with new values 
	delete[] state;
	state = NULL;

		
}

void Distribute2(int irow, int jcol,  double num)
{//distribute the Bugs according to beta with parameter gamma and nu
	if (num == 0) //if there are no Bugs don't worry
	{
		return;
	}
	else
	{
		double P_Sum = 0; // to keep track of probability used so that we readjust
		int numSum = 0;
		std::vector<double> MyDist;
		//copy the distribution vector to MyDist
		FlightP::GetBugFlight(MyDist);

		int NRows = Grid::GetNumRows();
		int NCols = Grid::GetNumCols();
		int Len = MyDist.size();


		//central point 

		double NewV(0);
		double myCrop = Grid::GetTreeDensity(irow, jcol); //get treedensity

		int N = 1;
		int M = num;
		double P = MyDist[0];
		int const LR = 1;
		double R[LR];
		int X[1];
		int ifail = 1;
		int MODE = 3;
		int const genid(1);
		int const subid(0);
		int lstate(0);
		int *state;
		state = new int[1];
		G05KGF(genid, subid, state, lstate, ifail);
		delete[] state;
		state = new int[lstate];
		ifail = 1;
		G05KGF(genid, subid, state, lstate, ifail);
		G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
		P_Sum = P_Sum + MyDist[0];
		M = M - X[0];


		if (myCrop>0.0) //Trees to land on
		{

			NewV = Grid::GetTempBug(irow, jcol) + X[0];
			Grid::SetTempBug(irow, jcol,  NewV);
			numSum = numSum + X[0];
		}
		else
		{
			flyAgain(irow, jcol);
			NewV = Grid::GetTempBug(irow, jcol) + X[0];
			Grid::SetTempBug(irow, jcol,  NewV);
			numSum = numSum + X[0];
		}

		int sideCount(2);
		int LenCount(3);
		while (Len>LenCount - 1) //next square%
		{
			//sides
			int nrow = irow;
			int ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);

			double myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop>0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
			}



			nrow = irow;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow + sideCount - 1;
			ncol = jcol;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //Get number in removed section

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Then there are trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - sideCount]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - sideCount];
			M = M - X[0];

			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
			}

			//corners

			nrow = irow + sideCount - 1;
			ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol);//get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0]; //alice check
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
			}
			//NewV=Grid::GetTempBug(nrow, ncol, health)+MyDist[LenCount-1]*num;
			//Grid::SetTempBug(nrow, ncol, health, NewV);

			nrow = irow + sideCount - 1;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'
			if (X[0] > 0)
			{
				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop > 0.0)//Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol + sideCount - 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol, NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
			}

			nrow = irow - sideCount + 1;
			ncol = jcol - sideCount + 1;
			Reflect(nrow, ncol);
			myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

			if (P_Sum < 1)
			{
				P = (MyDist[LenCount - 1]) / (1 - P_Sum);
			}
			else
			{
				P = 1;
			}
			if (P>1)
				P = 1;
			G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
			P_Sum = P_Sum + MyDist[LenCount - 1];
			M = M - X[0];
			if (X[0] > 0)
			{
				if (myCrop > 0.0) //Trees to land on
				{
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
				else
				{
					flyAgain(nrow, ncol);
					NewV = Grid::GetTempBug(nrow, ncol) + X[0];
					Grid::SetTempBug(nrow, ncol,  NewV);
					numSum = numSum + X[0];
				}
			}
			//diagonalsides

			int Mid = sideCount - 2;

			for (int Mcount = 0; Mcount<Mid; Mcount++)
			{
				nrow = irow + Mcount + 1;
				ncol = jcol + sideCount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}

				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol, NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - Mcount - 1;
				ncol = jcol + sideCount - 1;
				Reflect(nrow, ncol);
				double myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol, NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow + Mcount + 1;
				ncol = jcol - sideCount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - Mcount - 1;
				ncol = jcol - sideCount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				////
				nrow = irow + sideCount - 1;
				ncol = jcol + Mcount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow + sideCount - 1;
				ncol = jcol - Mcount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //Get proportion removed

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				//NewV=Grid::GetTempBug(nrow, ncol, health)+MyDist[LenCount-sideCount+Mcount+1]*num;
				//Grid::SetTempBug(nrow, ncol, health, NewV);
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //There is crop to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - sideCount + 1;
				ncol = jcol + Mcount + 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;

				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
				}

				nrow = irow - sideCount + 1;
				ncol = jcol - Mcount - 1;
				Reflect(nrow, ncol);
				myCrop = Grid::GetTreeDensity(nrow, ncol); //get proportion 'removed'

				if (P_Sum < 1)
				{
					P = (MyDist[LenCount - sideCount + Mcount + 1]) / (1 - P_Sum);
				}
				else
				{
					P = 1;
				}
				if (P>1)
					P = 1;
				G05TAF(MODE, N, M, P, R, LR, state, X, ifail);
				P_Sum = P_Sum + MyDist[LenCount - sideCount + Mcount + 1];
				M = M - X[0];
				if (X[0] > 0)
				{
					if (myCrop > 0.0) //Trees to land on
					{
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol,  NewV);
						numSum = numSum + X[0];
					}
					else
					{
						flyAgain(nrow, ncol);
						NewV = Grid::GetTempBug(nrow, ncol) + X[0];
						Grid::SetTempBug(nrow, ncol, NewV);
						numSum = numSum + X[0];
					}
				}


			}


			sideCount++; //these keep track of where you are on square see Fig 3 ECBSpatialModel.pdf
			LenCount = LenCount + sideCount;

		}
		delete[] state;
		state = NULL;
		if (numSum < num)
		{
			NewV = Grid::GetTempBug(irow, jcol) + num - numSum;
			Grid::SetTempBug(irow, jcol,  NewV);
		}
	}

		

		
}

void Reflect(int& nrow, int& ncol)
{
	int NRows=Grid::GetNumRows();
	int NCols=Grid::GetNumCols(); 
	while((nrow>NRows-1)||(nrow<0))
	{
		if (nrow<0)
		{
		nrow=-1-nrow;
		}
		else if (nrow>NRows-1)
		{
		nrow=2*NRows-1-nrow;
		}
	}
	while((ncol>NCols-1)||(ncol<0))
	{
		if (ncol<0)
		{
			ncol=-1-ncol;
		}
		else if
		(ncol>NCols-1)
		{
		ncol=2*NCols-1-ncol;
		}
	}


}

void flyAgain(int& nrow, int& ncol)
{
	//Find a near by field of maize
	
	int maxDis=Grid::GetNumCols();
	if (maxDis<Grid::GetNumRows()) maxDis=Grid::GetNumRows();
	for (int kcount=1; kcount<maxDis; kcount++)
	for (int irow=nrow-kcount; irow<nrow+kcount+1; irow++)
	{
		for (int jcol=ncol-kcount; jcol<ncol+kcount+1; jcol++)
		{
			if ((irow<Grid::GetNumRows())&&(irow>-1)&&(jcol<Grid::GetNumCols())&&(jcol>-1))
			{
				double myCrop = Grid::GetTreeDensity(irow, jcol);
				//if (Grid::GetCrop(irow, jcol, Grid::GetMaxNumStatus()-1)<0.999999)
				if (myCrop >0.0)
				{
					nrow=irow;
					ncol=jcol;
					return;
				}
			}

		}
	}
	
	throw std::logic_error("No Oranges for the Bugs :-("); 


}

void WriteCrops(int iyr)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char Crop[50];
    int n=sprintf_s(Crop, 50, "\\OutFiles\\TreeH%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutF(Npath_buffer);


	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeE%d.txt", iyr); //infected not infectious
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFE(Npath_buffer);
	
	
	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeI%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFI(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	n=sprintf_s(Crop, 50, "\\OutFiles\\TreeS%d.txt", iyr);
	strcat_s(Npath_buffer, Crop);
	std::ofstream OutFS(Npath_buffer);



	if ((OutF)&&(OutFI)&&(OutFS))
	{
		
			for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
			{
				for (int icount=0; icount<Grid::GetNumCols(); icount++)
				{
					OutF<<Grid::GetCrop(jcount, icount, 0)<<'\t';
					OutFE<<Grid::GetCrop(jcount, icount, 1)<<'\t';
					OutFI<<Grid::GetCrop(jcount, icount, 2)<<'\t';
					OutFS<<Grid::GetCrop(jcount, icount, 3)<<'\t';

				}
				OutF<<'\n';
				OutFE<<'\n';
				OutFI<<'\n';
				OutFS<<'\n';
			}
			
		
	}
	else
		throw std::logic_error("Error in write crops");

}

void SaveDataInStruc(double SaveStructure[], double SaveStructureT[], int tcount)
{
	int Nrows = Grid::GetNumRows();
	int Ncols = Grid::GetNumCols();
	int countcell = 0;
	for (int irow = 0; irow < Nrows; irow++)
	{
		for (int jcol = 0; jcol < Ncols; jcol++)
		{
			double myTreeDensity = Grid::GetTreeDensity(irow, jcol);
			double myIniTree = Grid::GetInitialTreeDensity(irow, jcol);
			if (myTreeDensity > 0)
			{
				// This needs to be some function of probability to detect	
				double tempTreeStates[12];
				for (int icount = 0; icount < Grid::GetMaxNumStatus(); icount++)
				{
					tempTreeStates[icount] = Grid::GetCrop(irow, jcol, icount); // These are proportions
				}

				//Get Infected trees 
				double propENoEAB = tempTreeStates[0] + tempTreeStates[1] + tempTreeStates[2] + tempTreeStates[3];
				if (propENoEAB > 1)
					propENoEAB = 1.0; // stops rounding errors 
				double propEInfGTYr = tempTreeStates[8] + tempTreeStates[9] + tempTreeStates[10] + tempTreeStates[11]; // GT a year of EAB
				double NumTreesHealthy = Grid::GetTreeDensity(irow, jcol) * propENoEAB;  // NB Tree Density is actually number !
				double NumTreesInfested = Grid::GetTreeDensity(irow, jcol) * (1 - propENoEAB);

				double myNumLarvae(0);
				double c_parameter(0.026); // Girdelled 0.5096; Traps 0.052; //visual = 0.026
				double NumTreeInspect(10);
				Grid::GetOWLav(irow, jcol, myNumLarvae);
				double phloemPerTree = Grid::GetPhloemPerTree();
				double myDetect = 0;

				/////////////////////////////////////////////////////////////////////////
				double P_LarvToAdult(0.0);
				Grid::GetPLarv(irow, jcol, P_LarvToAdult);
				///////////////////////////////////////////////////////////////////////
				double AveLavPerTree = myNumLarvae * phloemPerTree * NumTreesInfested / (myTreeDensity * P_LarvToAdult);
				int indx = tcount * Nrows * Ncols + irow * Ncols + jcol;
				SaveStructure[indx] = AveLavPerTree;
				SaveStructureT[indx] = myTreeDensity;
			}
			else
			{
				int indx = tcount * Nrows * Ncols + irow * Ncols + jcol;
				SaveStructure[indx] = 0;
				SaveStructureT[indx] = 0;
			}
		}
	}

}

void WriteDetectionInforOpt_Track(int occassion, int infRow, int infCol, double myWeight, int MaxYrs, double SaveStructure[], double SaveStructureT[])
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	strcpy_s(ECBStr, "\\OutFiles\\DataForOpt.txt");
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutF(Npath_buffer, std::ifstream::app);
	if (OutF)
	{
		int Nrows = Grid::GetNumRows();
		int Ncols = Grid::GetNumCols();


		for (int irow = 0; irow < Nrows; irow++)
		{
			for (int jcol = 0; jcol < Ncols; jcol++)
			{

				double myIniTree = Grid::GetInitialTreeDensity(irow, jcol);
				if (myIniTree > 0)
				{
					double scount = .0;
					for (int iicount = 0; iicount < MaxYrs; iicount++)
					{
						int indx = iicount * Nrows * Ncols + irow * Ncols + jcol;
						scount = scount + SaveStructure[indx];
					}
					if (scount > 0)
					{
						double myIniTree = Grid::GetInitialTreeDensity(irow, jcol);
						//OutF << infRow * Ncols + infCol << '\t' << myWeight << '\t' << irow * Ncols + jcol;
						OutF << occassion <<'\t' << infRow << '\t' << infCol << '\t' << myWeight << '\t' << myIniTree << '\t' << irow << '\t' << jcol;

						for (int tcount = 0; tcount < MaxYrs; tcount++)
						{
							int indx = tcount * Nrows * Ncols + irow * Ncols + jcol;
							OutF << '\t' << SaveStructure[indx];
						}
						for (int tcount = 0; tcount < MaxYrs; tcount++)
						{
							int indx = tcount * Nrows * Ncols + irow * Ncols + jcol;
							OutF << '\t' << SaveStructureT[indx];
						}
						OutF << '\n';
					}
					//else
					//{
						//throw std::logic_error("Bum2");
					//}
				}
				//else
				//{
					//throw std::logic_error("Bum");
				//}

			}
		}
		OutF.flush();
		OutF.close();
	}


}

void WriteLarvT(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	int n = sprintf_s(ECBStr, 50, "\\OutFiles\\LarvsT%d.txt", imth);
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutF(Npath_buffer);
	if (OutF)
	{
		for (int jcount = 0; jcount < Grid::GetNumRows(); jcount++)
		{
			for (int icount = 0; icount < Grid::GetNumCols(); icount++)
			{

				double totPop(0);
				 Grid::GetOWLav(jcount, icount, totPop ); //healthy 
				OutF << totPop << '\t';
				OutF.flush();

			}
			OutF << '\n';
		}
	}
	else
		throw std::logic_error("Error in write Larvae");

}

void WriteBugsT(int imth)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd( path_buffer, _MAX_PATH );
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
    int n=sprintf_s(ECBStr, 50, "\\OutFiles\\BugsT%d.txt", imth);
	strcat_s(Npath_buffer, ECBStr);
		
	std::ofstream OutF(Npath_buffer);
	if (OutF)
	{
		for (int jcount=0; jcount<Grid::GetNumRows(); jcount++)
		{
			for (int icount=0; icount<Grid::GetNumCols(); icount++)
			{
		
			    double totPop;
				totPop=Grid::GetNumG(jcount, icount ); //healthy 
				OutF<<totPop<<'\t';
				OutF.flush();

			}
			OutF<<'\n';
		}
	}
	else
		throw std::logic_error("Error in write bugs");

}


double ProbabilityofDetection(int LookHere[][ConstDefNumLocs], int NumLocs)
{
	double ProbDetect(0);

	for (int icount = 0; icount < NumLocs; icount++)
	{
		ProbDetect= ProbDetect + Grid::GetProbDetect(LookHere[0][icount], LookHere[1][icount]); //Number of Bugs
		
		//Grid::GetOWLav(LookHere[0][icount], LookHere[1][icount], totPop); //healthy 

	}

	return ProbDetect;


}


void DetectionFullRunW_Target(double RandIf, double kill_rate,  int ConstMaxInvYr)
{
	int iflag(0);
	//This function reads in the entry locations and runs the simulation for 8 years
	//The output is used in probability of detection in simulated annealing in a seperate code
	
	int MaxGridS = Grid::GetNumRows() * Grid::GetNumCols() * ConstMaxInvYr;
	double* mySaveStructure;
	double* mySaveStructureTree;
	mySaveStructure = new double[MaxGridS];
	mySaveStructureTree = new double[MaxGridS];

	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	strcpy_s(ECBStr, "\\TreeDensity\\SampleEntry.txt"); //This is a list of likely entry points. They resampled so repeat
	strcat_s(Npath_buffer, ECBStr);
	int myX, myY;
	int Occassion(0);
	
	std::ifstream InFme(Npath_buffer);
	if(InFme)
	{ 
		int LazyAlice = 10000; //length of sample input file
		for (int jcount = 0; jcount < LazyAlice; jcount++)
		{
				InFme>>myX;
				InFme>>myY;

				double TempTreeDensity = Grid::GetInitialTreeDensity(myY, myX);
				if (TempTreeDensity > 0)
				{
					Grid::ReSetLandBasic(myY, myX, iflag);

					double myWeight = Grid::GetInclustionProb(myY, myX);
					if (myWeight > 0) // if there is a chance this could be an entry point then run - else don't bother This is now a bug check!
					{
						for (int tcount = 0; tcount < ConstMaxInvYr; tcount++)
						{
							ModelMonth(tcount, RandIf, kill_rate);
							// create function to set probability of detection for each cell year
							SaveDataInStruc(mySaveStructure, mySaveStructureTree, tcount);
						}
						WriteDetectionInforOpt_Track(Occassion, myY, myX, myWeight, ConstMaxInvYr, mySaveStructure, mySaveStructureTree);
						Occassion = Occassion + 1;
					}
					else
						throw std::logic_error("Risk is zero please check");
				}
			
		}


	}
	else
		throw std::logic_error("Error in Read Samples");

	delete[] mySaveStructure;
	mySaveStructure = NULL;
	delete[] mySaveStructureTree;
	mySaveStructureTree = NULL;



}

void WriteAveLarvae(int imth, int Repeat)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	int n = sprintf_s(ECBStr, 50, "\\OutFiles\\AveLarvae%d_%d.txt",Repeat, imth);
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutF(Npath_buffer);

	double phloemPerTree = Grid::GetPhloemPerTree();
	if (OutF)
	{
		for (int icount = 0; icount < Grid::GetNumRows(); icount++)
		{
			for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
			{		
				double AveLavPerTree = 0.0;
				
				if (Grid::GetTreeDensity(icount, jcount) > 0.0)
				{
					double TreeDensity = Grid::GetTreeDensity(icount, jcount);
					// infected
					double tempTreeStates[12];
					for (int ii = 0; ii < Grid::GetMaxNumStatus(); ii++)
					{
						tempTreeStates[ii] = Grid::GetCrop(icount, jcount, ii); // These are proportions
					}

					//Get Infected trees 
					double propENoEAB = tempTreeStates[0] + tempTreeStates[1] + tempTreeStates[2] + tempTreeStates[3];
					if (propENoEAB > 1)
						propENoEAB = 1;
					double NumTreesInfested = TreeDensity * (1.0 - propENoEAB);

					double myNumLarvae(0.0);
					Grid::GetOWLav(icount, jcount, myNumLarvae);

					double P_LarvToAdult(0.0);
					Grid::GetPLarv(icount, jcount, P_LarvToAdult);
					AveLavPerTree = myNumLarvae * phloemPerTree * NumTreesInfested / (TreeDensity * P_LarvToAdult);
				}

					 
				OutF << AveLavPerTree << '\t';
				OutF.flush();

			}
			OutF << '\n';
		}
	}
	else
		throw std::logic_error("Error in write AveLarvaePerTree");

}

void WritePropDead(int imth, int Repeat)
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	int n = sprintf_s(ECBStr, 50, "\\OutFiles\\PropDeadAsh%d_%d.txt",Repeat, imth);
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutF(Npath_buffer);

	if (OutF)
	{
		for (int icount = 0; icount < Grid::GetNumRows(); icount++)
		{
			for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
			{
				double PropDead = 0.0;

				if (Grid::GetInitialTreeDensity(icount, jcount) > 0.0)
				{
					double TreeDensity = Grid::GetTreeDensity(icount, jcount);
					double InitialDensity = Grid::GetInitialTreeDensity(icount, jcount);

					PropDead = 1- TreeDensity / InitialDensity;
				}


				OutF << PropDead << '\t';
				OutF.flush();

			}
			OutF << '\n';
		}
	}
	else
		throw std::logic_error("Error in write ProportionAshDead");

}

void WriteProbDetectForR()
{
	char path_buffer[_MAX_PATH];
	char Npath_buffer[_MAX_PATH];
	_getcwd(path_buffer, _MAX_PATH);
	strcpy_s(Npath_buffer, path_buffer);
	char ECBStr[50];
	strcpy_s(ECBStr, 50, "\\OutFiles\\ProbDetectForR.txt");
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutF(Npath_buffer);

	strcpy_s(Npath_buffer, path_buffer);
	strcpy_s(ECBStr, 50, "\\OutFiles\\ProbDetectForMat.txt");
	strcat_s(Npath_buffer, ECBStr);

	std::ofstream OutFF(Npath_buffer);

	if (OutF)
	{
		for (int icount = 0; icount < Grid::GetNumRows(); icount++)
		{
			for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
			{
				double myprobDetect=Grid::GetProbDetect(icount, jcount);

				OutF <<icount << '\t' << jcount << '\t'<< myprobDetect<<'\n';
				OutF.flush();
				OutFF <<  myprobDetect << '\t';
				OutFF.flush();

			}
			OutFF << '\n';
		}
	}
	else
		throw std::logic_error("Error in write ProbDetectForR");

}