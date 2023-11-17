#ifndef __EABGrid_H__
#define __EABGrid_H__

#include <vector>
#include "D:\\Users\\milne\\AppData\\Local\\NAG\\NL29\\nlw3229del\\include\\nag.h"

namespace Grid
{
	int GetNumRows();
	int GetNumCols();
	double GetPropOverWinterLarvaeStart();

	void SetEmergenceDensity(int irow, int jcol, double Stype);
	double GetEmergenceDensity(int irow, int jcol);

	//Getting and setting numbers of HLB
	double GetNumG(int irow, int jcol); //get total number of adults in a cell
	void SetNumG(int irow, int jcol, double NewNum); //set total number of bugs
	
	void SetTempBug(int irow, int jcol,  double NewNum);
	double GetTempBug(int irow, int jcol);
	double GetGridLen();
	
	//Keeps a vector of larvae MAY NOT NEED
	void SetOWLav(int irow, int jcol, double totPopp); //get numbers of Bugs after they have grown and moved
	void GetOWLav(int irow, int jcol,  double& totPop); //get numbers of Bugs

	void GetPLarv(int irow, int jcol, double& P_LarvToAdult); // Proportion of OW larvae 
	void SetPLarv(int irow, int jcol, double P_LarvToAdult); // Proportion of OW larvae
	
	//Seting trees status
	void SetCrop(int irow, int jcol, int myStat, double per); // sets proportion of trees at a given status
	double GetCrop(int irow, int jcol, int myStat);
	int GetMaxNumStatus(); // retuns the number of status a tree can be in , i.e. health, infected, sympomatic, removed
	
	//These are related to number of cells in plantations might need revising for our needs
	int GetNumCells(int PlantaCount); //returns number of grid cells in a plantation. Probably always 4
	int GetCellLocation(int PlantaCount, int cellCount); // thsi returns a number representing location icol*numRows+irow
	
	//Set up landscape comands
	void SetLand(int newL, int& iflag); 
	void SetLandBasic(int& iflag);
	void ReSetLandBasic(int infCelljRow, int infCelliCol, int& iflag);
		
	void SetPlanta(int irow, int jcol, int Ctype);
	void SetTreeDensity(int irow, int jcol, double Ctype); // This is now tree density needs renaming SetCHMA
	void SetInitialTreeDensity(int irow, int jcol, double Ctype); // This is now tree density needs renaming SetCHMA
	//Get landscape info
	
	double GetTreeDensity(int irow, int jcol);  // This is now tree density needsGetCHMA
	double GetInitialTreeDensity(int irow, int jcol);  // This is now tree density needsGetCHMA
	
	
	//Apply felling or not
	void GetTreeFellingStats(int icount, int jcount, int& FellingType, double& propFelled); // Felling type:  0 = no felling; 1 = fell EAB infected; 2= sanitation felling; propFelled tells the proportion of ill trees to fell
	void SetTreeFellingStats(int icount, int jcount, int FellingType, double propFelled); // Felling type:  0 = no felling; 1 = fell EAB infected; 2= sanitation felling; propFelled tells the proportion of ill trees to fell

	void SetInclustionProb(int jrow, int icol, double IncPr); //Probability that this is the entry location of the first infection 
	double GetInclustionProb(int jrow, int icol); //Probability that this is the entry location of the first infection 

	void ClearData();
	void ClearTempBugs();
	
	void SetIntVars(double jdist, double gridLen, double lambda);
	void GetIntVars(double& jdist, double& gridLen, double& lambda);
	
		
	//Comands to track invasion may not need
	void SetInvasion(int irow, int jcol, int yr);
	int GetInvasion(int irow, int jcol);
	double GetPhloemPerTree();
	void SetTotalPlantasFilled(int mycount);
	int GetTotalPlantasFilled();

	//OptimisationStff
	void ReSetProbDetection();
	void SetProbDetection();
	void SetWeightedProbDetection(double myWeight); // DO NOT MIX 
	double GetProbDetect(int irow, int jcol);
	
}

namespace FlightP
{
	void ClearData();
	void Initialise();
	void GetBugFlight(std::vector<double>& MyDist);


}

double Integrate(int idist, int jdist, double lambda);


void GridIDToCoord(int Floc, int& irow, int& icol);
int CoordToGridID(int irow, int icol);

//void 
extern "C" void __stdcall fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] );


#endif