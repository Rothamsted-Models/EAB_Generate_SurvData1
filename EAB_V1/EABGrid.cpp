#include "EABGrid.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <direct.h>


//int const NumCropTypes=4;
int const MaxNumNeighbours=50;
int const MaxNumPCells=12169; //maximum number cells per plantation



//typedef double ResHealthSex[2];// index: 0=healthy, 1=infected - not needed
	
typedef double ThreeDoubles[3];

class PInfo
{
public:
	PInfo()
	{
	//	BlockNum=-1;
		DisStatus=0; //We assume the worry level is zero on initialisation
		Score=0;
		BeliefInControl = 0; 
		MonthlyLossofHealth = 0;
	//	numCells=0;
		//for (int icount=0; icount<MaxNumPCells; icount++)
		//{
		Cells.clear();
		//}

		for (int icount=0; icount<MaxNumNeighbours; icount++)
		{
			NeighbourBlockID[icount]=-1; // list of neighbours
		}
		NumNeighbours=0; //number of neighbours
		ConnectionType=-1;
		UniquePlantaID=-1;
		

	}
	PInfo(int myPlantaID,  int mySym, double myscore, double mybelief, std::vector<int> mycells)
	{
		//FCHMA=myCHMA; // now tree density needs renaming
	//	BlockNum=-1; //not used
		UniquePlantaID=myPlantaID;
		DisStatus=mySym; //0 if < 5% sympotmatic and 1 if more than that in any one cell
		Score=myscore;
		BeliefInControl = mybelief;
		MonthlyLossofHealth = 0;
		//numCells=myNumcells;
		Cells.clear();
		for (int icount=0; icount<mycells.size(); icount++)
		{
			Cells.push_back(mycells[icount]);
		}
		for (int icount=0; icount<MaxNumNeighbours; icount++)
		{
			NeighbourBlockID[icount]=-1; // list of neighbours
		}
		NumNeighbours=0; //number of neighbours
		ConnectionType=-1;
		
	}
	~PInfo()
	{
		NumNeighbours=0; //number of neighbours
		Cells.clear();
	}
	
	void SetPBelief(double mybelief)
	{
		BeliefInControl = mybelief;
	}
	void SetNeighbour(int NID)
	{
		if (NumNeighbours<MaxNumNeighbours) 
		{
			bool alreadyCounted=false;
			for (int icount=0; icount<NumNeighbours; icount++)
			{
				if (NeighbourBlockID[icount]==NID)
						alreadyCounted=true;
			}
			if (!alreadyCounted)
			{
				NeighbourBlockID[NumNeighbours]=NID;
				NumNeighbours++;
			}
		}
		else
		{
			throw std::logic_error("too many neighbours");
		}
	}
	void SetConnect(int myCon)//This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
	{
		ConnectionType=myCon;
	}
	
	int GetnumCells(){return Cells.size();}
	int GetCellID(int CellCount){return Cells[CellCount];}
	int GetnumNeighbours(){return NumNeighbours;}
	int GetNeighbourID(int ncount){return NeighbourBlockID[ncount];}
	double GetScore(){return Score;}
	double GetBelief(){ return BeliefInControl; }
	//double GetFmCHMA(){return FCHMA;}
//	int GetFmPlanta(){return UniquePlantaID;}
	void GetNeighbour(int NID, int myNeighbs[], int& NumNeibs)
	{
		for (int icount=0; icount<NumNeighbours; icount++)
		{
			myNeighbs[icount]=NeighbourBlockID[icount];
		}

		NumNeibs=NumNeighbours;
		
	}
	int GetPlantNum(){return UniquePlantaID;}
	void DeleteNeigh(){NumNeighbours=0;}
	int GetConnect()//This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
	{
		return ConnectionType;
	}
	int GetDisStatus(){ return DisStatus;}
	void SetDisStatus(int myStat){DisStatus=myStat;}
	void AddHLoss(double exloss){ MonthlyLossofHealth = MonthlyLossofHealth + exloss; }
	double GetHLoss(){ return MonthlyLossofHealth; }
	

private:
	//	int BlockNum;
		int DisStatus; //If symptoms above 5% somewhere on the plantation then worry sets in
		double Score; //this is the 'worry' 0 is not worried about infection and 1 is very worried
		double BeliefInControl; //If a farmer firmly believes in the control strategy then we score him 1 is not zero. 
		std::vector<int> Cells; // field locations field=icols*maxnumrows+irows
		//int numCells; // number of cells in a Block
		//double  FCHMA; 
		int NeighbourBlockID[MaxNumNeighbours]; // list of neighbours
		int NumNeighbours; //number of neighbours
		int ConnectionType; //This tells us how the Block is connected to the information highway. -1 = no connection, 0=neighbour; 1=Planta; 2=CHMA
        int UniquePlantaID; //the Planta the Block is in
		double MonthlyLossofHealth;
		
};

class TwoVars
{
public:
  TwoVars(){Xs[0]=0; Xs[1]=0;};
  TwoVars(double val1, double val2){Xs[0]=val1; Xs[1]=val2;};
  ~TwoVars(){};
  TwoVars(const TwoVars& ATheX){Xs[0]=ATheX.Xs[0]; Xs[1]=ATheX.Xs[1];};
  TwoVars& operator=(const  TwoVars& ATheX){Xs[0]=ATheX.Xs[0]; Xs[1]=ATheX.Xs[1]; return *this;};
  bool  operator<(const  TwoVars& ATheX)
  {
    bool A=false; 
    if(Xs[0]<ATheX.Xs[0]) A=true;
    return A;
  };
  bool  operator>(const  TwoVars& ATheX){bool A=false; if(Xs[0]>ATheX.Xs[0]) A=true; return A;};
  double Xs[2];
};

class Bugs
{
public:
	
	Bugs()
	{
			BugT=0;
	}

	~Bugs()
	{
		TotalOWLav = 0;
		BugT=0;
	}
	//double GetNum(int Health){return BugT[Health];} //get total number of specified type
	double GetNum() { return BugT; } //get total number of specified type
		
	//void SetNum(int Health, double NewNum){BugT[Health]=NewNum;} 
	void SetNum( double NewNum) { BugT = NewNum; }

	void SetOWLav(double totPop) { TotalOWLav = totPop; }
	
	void GetOWLav(double& totPop){totPop=TotalOWLav;}
	
	void ClearTempBug()
	{
		//TempBug[0]=0; //healthy
		//TempBug[1]=0; //sick
		TempBug=0;
		
	}
	void Clear()
	{
		
		TotalOWLav=0;
		BugT = 0;
		TempBug = 0;

	}
	void SetTempBug(  double NewNum)
	{
		TempBug=NewNum;
	}
	double GetTempBug()
	{
		return TempBug;
	}


private:
	double BugT; // number of each type of adult Bug (average per plant/unit area)
	double TotalOWLav;
	double TempBug; //a structure to keep track of travelling Bugs
	
};



namespace Grid
{
	int const HistLen(12); // This is how far back growers look when assessing evidence. We may not need it
	//int const MaxNumYears(200);

	
	
	//Size of grid - we should make a 100 x 100 one to start
	int const Nrows = 1212;   //UK 1212;//SE 275; // 100; //<-Test // 515; //<-UKGrid //515; //Test // 100; // dummy //265;  //CHMA 54 //
	int const Ncols = 602;   //UK 602;//SE206; // 100; //<-Test //535; //<-UKGrid //535; //Test // 100; // dummy //220; //CHMA 54 //
	int const NStatus=12; //Status of trees. 0-3 not infested by EAB; 4-7 infested last year; 8-11 infested before last year
	double const gridLen = 1.0; // grid length is 1 = 1km
	double const GridVal_P_LarvToAdultStart = 0.25; //0.95;//LarvaeToAdultInOneYear(Lyons, 2015) // set tp 0.95
	double const TreesPerHa = 889; // Trees per ha informed from the National Forest Inventory epimath:\EAB\Project Documents\Population Dynamics
	double const G_phloemPerTree = 4.53765; //m2 of pholem per tree - very estimated needs checking - if you wabt to do as pholem change to 1 epimath:\EAB\Project Documents\Population Dynamics
	double GetGridLen() { return gridLen; };
	double GP_LarvToAdult;


	

	// This is either Healthy=0, 17
	//No EAB +Exposed  by ADB but not infectious = 1, 
	//No EAB +infected by ADB infectctious no symtoms=2, 
	//No EAB + Symptoms by ADB =3, 
	//EAB but not detected + no ADB =4, 
	//EAB but not deteted +Exposed by ADB but not infectious = 5, 
	//EAB but not deteted +infected by ADB infectctious no symtoms=6, 
	// EAB but not deteted + Symptoms by ADB =7, 
	//EAB detected + no ADB =8, 
	//EAB deteted +Exposed by ADB but not infectious = 9, 
	//EAB  deteted +infected by ADB infectctious no symtoms=10, 
	// EAB deteted + Symptoms by ADB =11, 
	double Trees[Nrows][Ncols][NStatus]; //Holds proportion of trees in each CHMA listed under status
	Bugs GridBugs[Nrows][Ncols];

	double ProbDetection[Nrows][Ncols]; // This holds the probability of detection for each location for each realisation and year
	
	double ProbInitialInfectionDetection[Nrows][Ncols]; // This holds the probability this site is the entry location 


	//May not be needed
	int G_Control[Nrows][Ncols][HistLen]; //We now hold control for 12 months. the zero position is control for last month the 1 position for the month before
	//Note that  areas that cannot be cropped are 4 with fix=1
	
	int G_FellType[Nrows][Ncols]; // This is the type of felling in the cell 0=none; 1 = a proportion of infested with EAB; 2 = sanitation felling
	double G_PropFelled[Nrows][Ncols]; // If G_FellType = 1 this tells you the prportion to fell


	double IntVars[3]; //keeps parameters for dispersal integration can be made bigger
	//each cell has a CHMA, Planta and Block ID. The combination is unique.
	double G_TreeDensity[Nrows][Ncols]; // HoldsTreeDensity so implicitly the number removed (trees per km square)
	double IniG_TreeDensity[Nrows][Ncols]; // Holds Initial TreeDensity  (trees per km square)
	double Emerge_m2[Nrows][Ncols]; // store values before flight!


	int PlantaID[Nrows][Ncols]; // Holds ownership might need to be more complex as likely to have multiple ownders per 1km cell
	//int BlockID[Nrows][Ncols];

	int invasionYr[Nrows][Ncols]; // this is to store the year that the invasion of the Bug becomes > a fixed number
	std::vector<PInfo> PlantationInfo; //This stores the info about each Block. Uniquie Block ID is offset by one from indexing, i.e. Block 1 info is held in index 0
	int TotalPlantasFilled;

	int GetTotalPlantasFilled(){return TotalPlantasFilled;}

	int GetNumRows(){return Nrows;}
	
	int GetNumCols(){return Ncols;}

	double GetPropOverWinterLarvaeStart() { return GridVal_P_LarvToAdultStart; };

	double GetPhloemPerTree() { return G_phloemPerTree; }

	void SetInclustionProb(int irow, int jcol, double IncPr) //Probability that this is the entry location of the first infection 
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			ProbInitialInfectionDetection[irow][jcol] = IncPr;
		}
	}
	double GetInclustionProb(int irow, int jcol) //Probability that this is the entry location of the first infection 
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			return ProbInitialInfectionDetection[irow][jcol];
		}
	}

	//int GetNumResTypes(){return NumResTypes;}

	void SetEmergenceDensity(int irow, int jcol, double Stype) // was SetCHMA
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			Emerge_m2[irow][jcol] = Stype;
		}
		else
			throw std::logic_error("Error in SetEmergeDensity");
		//int junk=1;

	}
	
	double GetEmergenceDensity(int irow, int jcol) // Now tree density GetCHMA
	{
		if ((irow > -1) && (irow < Nrows) && (jcol > -1) && (jcol < Ncols))
			return Emerge_m2[irow][jcol];
		else
			throw std::logic_error("error in get Emergence Density");
		//int junk=-9;
	}

    double GetNumG(int irow, int jcol)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			return GridBugs[irow][jcol].GetNum();
		}
		else
			throw std::logic_error("Error in GetNumG");
	
	} //get total number of specified type
	
	
	int GetNumCells(int BlockID){return PlantationInfo[BlockID].GetnumCells();}

	int GetCellLocation(int PlantaID, int cellCount){ return PlantationInfo[PlantaID].GetCellID(cellCount);}

	void SetTotalPlantasFilled(int myCounts){ TotalPlantasFilled=myCounts;}

	void SetNumG(int irow, int jcol, double NewNum)
	{ 
		if ((irow<GetNumRows())&&(jcol<GetNumCols()))
		{
			GridBugs[irow][jcol].SetNum(NewNum);
		}
		else
			throw std::logic_error("Error initialising in grid cells that don't exist");

	} 

	
	void SetOWLav(int irow, int jcol, double totPop){ GridBugs[irow][jcol].SetOWLav(totPop);}
	
	
	

	void SetPlanta(int irow, int jcol, int Ctype)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			PlantaID[irow][jcol]=Ctype; 
		}
		else
			throw std::logic_error("Error in SetPlanta");
			//int junk=1;

	}


	void GetPLarv(int irow, int jcol, double& P_LarvToAdult)
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			P_LarvToAdult=GP_LarvToAdult ;
		}
		else
			throw std::logic_error("Error in GetPLarv");
	}
	void SetPLarv(int irow, int jcol, double P_LarvToAdult)
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			 GP_LarvToAdult= P_LarvToAdult;
		}
		else
			throw std::logic_error("Error in GetPLarv");
	}


	
	void SetTreeDensity(int irow, int jcol, double Stype) // was SetCHMA
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1))
		{
			G_TreeDensity[irow][jcol]=Stype;
		}
		else
			throw std::logic_error("Error in SetTreeDensity");
			//int junk=1;

	}

	void SetInitialTreeDensity(int irow, int jcol, double Stype) // was SetCHMA
	{
		if ((irow < GetNumRows()) && (jcol < GetNumCols()) && (irow > -1) && (jcol > -1))
		{
			IniG_TreeDensity[irow][jcol] = Stype;
		}
		else
			throw std::logic_error("Error in SetTreeDensity");
		//int junk=1;

	}
	
	void SetCrop(int irow, int jcol, int ist, double per)
	{
		if ((irow<GetNumRows())&&(jcol<GetNumCols())&&(irow>-1)&&(jcol>-1) &&(ist<NStatus))
		{

			if ((per < -0.00000001) || (per > 1.000000001)) { // allowing rounding error
			//	int junk = 1;
				throw std::logic_error("Error in SetCrop1");
			}
			if ((per > 1)) // if a little over pull it back
				per = 1.0;
			if ((per < 0))
				per = 0.0;
			Trees[irow][jcol][ist]=per; 
			
		}
		else
			throw std::logic_error("Error in SetCrop2");
			//int junk=1;

	}

	

		
	
	double GetCrop(int irow, int jcol, int ist)
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return Trees[irow][jcol][ist];
		else
			throw std::logic_error("error in get crop");
			//return -9;
	}

	

	double GetTreeDensity(int irow, int jcol) // Now tree density GetCHMA
	{
		if ((irow>-1)&&(irow<Nrows)&&(jcol>-1)&&(jcol<Ncols))
			return G_TreeDensity[irow][jcol];
		else
			throw std::logic_error("error in get TreeDensity");
			//int junk=-9;
	}

	double GetInitialTreeDensity(int irow, int jcol) // Now tree density GetCHMA
	{
		if ((irow > -1) && (irow < Nrows) && (jcol > -1) && (jcol < Ncols))
			return IniG_TreeDensity[irow][jcol];
		else
			throw std::logic_error("error in get TreeDensity");
		//int junk=-9;
	}


	int GetMaxNumStatus(){return NStatus;}
			
	void GetOWLav(int irow, int jcol,  double& totHPop){ GridBugs[irow][jcol].GetOWLav(totHPop);}

	
	
	
	void ClearData()
	{
		for (int icount=0; icount<Nrows; icount++)
		{
			for (int jcount=0; jcount<Ncols; jcount++)
			{
				GridBugs[icount][jcount].Clear();
				for (int scount=0; scount<NStatus; scount++)
				{
					Trees[icount][jcount][scount]=0;
				}
			}
		}

		PlantationInfo.clear();

	}

	void ClearTempBugs()
	{
		for (int icount=0; icount<Nrows; icount++)
		{
			for (int jcount=0; jcount<Ncols; jcount++)
			{
				GridBugs[icount][jcount].ClearTempBug();
				
			}
		}
	}

	void SetTempBug(int irow, int jcol,  double NewNum)
	{
		GridBugs[irow][jcol].SetTempBug( NewNum);

	}
	double GetTempBug(int irow, int jcol)
	{
		if ((irow>-1)&&(irow<Grid::GetNumRows())&&(jcol>-1)&&(jcol<Grid::GetNumCols()))
			return GridBugs[irow][jcol].GetTempBug();
		else
			throw std::logic_error("index out of bounds");
	}

	
	//Apply felling or not
	void GetTreeFellingStats(int irow, int jcol, int& FellingType, double& propFelled) // Felling type:  0 = no felling; 1 = fell EAB infected; 2= sanitation felling; propFelled tells the proportion of ill trees to fell
	{
		if ((irow > -1) && (irow < Grid::GetNumRows()) && (jcol > -1) && (jcol < Grid::GetNumCols()))
		{
			FellingType = G_FellType[irow][jcol];
			propFelled  = G_PropFelled[irow][jcol];
		}

		else
			throw std::logic_error("Error in GetTreeFellingStats index out of bounds");

	}
		
	void SetTreeFellingStats(int irow, int jcol, int FellingType, double propFelled) // Felling type:  0 = no felling; 1 = fell EAB infected; 2= sanitation felling; propFelled tells the proportion of ill trees to fell
	{
		if ((irow > -1) && (irow < Grid::GetNumRows()) && (jcol > -1) && (jcol < Grid::GetNumCols()))
		{
			G_FellType[irow][jcol]= FellingType;
			G_PropFelled[irow][jcol]= propFelled;
		}

		else
			throw std::logic_error("Error in SetTreeFellingStats index out of bounds");

	}

	void SetIntVars(double jdist, double gridLen, double lambda)
	{
		IntVars[0]=jdist;
		IntVars[1]=gridLen;
		IntVars[2]=lambda;

	}

	void GetIntVars(double& jdist, double& gridLen, double& lambda)
	{
		jdist=IntVars[0];
		gridLen=IntVars[1];
		lambda=IntVars[2];
		
	}

	
	void SetInvasion(int irow, int jcol, int yr) //stores where Bugs have reached in an invasion from left to right
	{
		invasionYr[irow][jcol]=yr;
	}

	int GetInvasion(int irow, int jcol)
	{
		return invasionYr[irow][jcol];
	}

	void SetLand(int newL, int& iflag) //double alpha_W, double beta_W, double alpha_B, double beta_B are paraters for beta dristribution describing initial belief structures 
	//void SetLand(int newL, double alpha_W, double beta_W, double alpha_B, double beta_B, int& iflag) //double alpha_W, double beta_W, double alpha_B, double beta_B are paraters for beta dristribution describing initial belief structures 
	{
		
		int UBlockID=1; //unique Block ID	

		for (int icount=0;icount<Nrows; icount++)
		{
			for (int jcount=0;jcount<Ncols; jcount++)
			{
				SetTreeDensity(icount,jcount,0.0); // This is just cleaning it before we start was SetCHMA(icount,jcount,0.0); 
				SetInitialTreeDensity(icount, jcount, 0.0); // This is just cleaning it before we start was SetCHMA(icount,jcount,0.0); 
				SetPlanta(icount,jcount,-9);
				for (int kcount=0; kcount<NStatus; kcount++)
					SetCrop(icount,jcount,kcount,0);  //set all status to zero
				
				
			}
		}
		if (newL==0) //then read grid spec from file 
		{
			char path_buffer[_MAX_PATH];
			char Npath_buffer[_MAX_PATH];
			_getcwd( path_buffer, _MAX_PATH );
			strcpy_s(Npath_buffer, path_buffer);
			//strcat_s(Npath_buffer, "\\TreeDensity\\UKGridOwner.txt"); 
			//strcat_s(Npath_buffer, "\\TreeDensity\\Owners_SE.txt");
			strcat_s(Npath_buffer, "\\TreeDensity\\Owners_GB.txt");

			std::ifstream    InF(Npath_buffer);

			strcpy_s(Npath_buffer, path_buffer);
			//strcat_s(Npath_buffer, "\\TreeDensity\\UKGridADB.txt"); 
			//strcat_s(Npath_buffer, "\\TreeDensity\\ADBInfection_SE.txt"); 
			strcat_s(Npath_buffer, "\\TreeDensity\\ADBInfection_GB.txt");
			std::ifstream    InFF(Npath_buffer);

			strcpy_s(Npath_buffer, path_buffer);
			//strcat_s(Npath_buffer, "\\TreeDensity\\UKGrid.txt"); 
			//strcat_s(Npath_buffer, "\\TreeDensity\\TreeDensity_SE.txt");
			strcat_s(Npath_buffer, "\\TreeDensity\\treeDensity_GB.txt");
			std::ifstream	 InFFF(Npath_buffer);

			strcpy_s(Npath_buffer, path_buffer);
			//Create files for SE; 
			
			//strcat_s(Npath_buffer, "\\TreeDensity\\InclusionProbSE.txt");
			strcat_s(Npath_buffer, "\\TreeDensity\\InclusionProbGB.txt");

			std::ifstream	 InF4(Npath_buffer);

			int TempOwner;
			double TempTreeDensity, ADBInf, infectionStatus;
			int Sick(0); //Flag to see of sick tree set/
			int iyr=2;
			int MaxPlantID=-1; // This keeps track of the maximum ID number for the plantations for use later 
			
			///if (InF && InFF && InFFF)
			if (InF&&InFF&&InFFF&&InF4)
			{
				//Sample from uniform
				int ifail = 1;
				int lstate(0);
				int subid(1);
				int* state;
				state = new int[lstate];
				int Lseed = 1;
				int seed[1];
				seed[0] = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail); //always get same vale of random permutaion for a given seed
				delete[] state;
				state = new int[lstate];
				ifail = 1;
				G05KGF(1, subid, state, lstate, ifail);
				//G05KFF(1,subid, seed, Lseed, state,lstate,ifail);
				double* urand;
				int NumVals = 1;
				urand = new double[NumVals];
				ifail = 1;
				G05SAF(NumVals, state, urand, ifail);
				if (ifail>0)
					throw std::logic_error("Error in sampling intial infection");
				double myCountProb = 0;
				int infCelliCol = 0;
				int infCelljRow = 0;
				urand[0] = urand[0] * 100.0; //turn to percents to matxh input file
				int notFoundIT = 1;

				for (int jcount = 0; jcount < Nrows; jcount++)
				{
					for (int icount = 0; icount < Ncols; icount++)
					{
						InF4 >> infectionStatus;
						myCountProb= myCountProb+ infectionStatus;
						if ((myCountProb> urand[0]) && (notFoundIT))
						{
							infCelliCol = icount;
							infCelljRow = jcount;
							notFoundIT = 0;
						}
					}
				}
				delete[]  urand;
				urand = NULL;
				delete[] state; 
				state = NULL;

				for (int jcount=0; jcount<Nrows; jcount++)
				{
					for (int icount=0; icount<Ncols; icount++)
					{
						InF>> TempOwner;
						PlantaID[jcount][icount]= TempOwner;
						if (TempOwner >MaxPlantID)
							MaxPlantID= TempOwner;
						InFF>>ADBInf; // This is assumed to be symtomatic we make no attempt to predict other states
						//BlockID[jcount][icount]=TempBlock;
						ADBInf = 0;
						InFFF>> TempTreeDensity;
						TempTreeDensity = TempTreeDensity*TreesPerHa; //This is number of trees per cell = ha / cell * trees / ha
						G_TreeDensity[jcount][icount]= TempTreeDensity;
						IniG_TreeDensity[jcount][icount] = TempTreeDensity;
						if (TempTreeDensity < 1.0)
							TempTreeDensity = 1;
						double ProbIniTreeInfection = 1 / TempTreeDensity; // This is now set so we have 1 tree infected per ha.
										
						//Here we set the first cell that has a crop to 10% EAB landed 
						if (TempTreeDensity >0)
						{
							SetCrop(jcount, icount, 0, 1 - ADBInf); //All healthy
							SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
							SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
							SetCrop(jcount, icount, 3, ADBInf); //ADB: symp
							SetCrop(jcount, icount, 4, 0.0); //EAB present not detected + No ADB
							SetCrop(jcount, icount, 5, 0.0); //EAB present not detected +ADB: Exposed not syp
							SetCrop(jcount, icount, 6, 0.0); //EAB present not detected +ADB: Infected 
							SetCrop(jcount, icount, 7, 0.0); //EAB present not detected +ADB: symp
							SetCrop(jcount, icount, 8, 0.0); //EAB detected + No ADB
							SetCrop(jcount, icount, 9, 0.0); //EAB detected +ADB: Exposed not syp
							SetCrop(jcount, icount, 10, 0.0); //EAB detected +ADB: Infected 
							SetCrop(jcount, icount, 11, 0.0); //EAB detected +ADB: symp

							if ((jcount == infCelljRow) && (icount == infCelliCol))
							{

								double EABonADB = (ProbIniTreeInfection)*ADBInf;
								double EABNotonADB = (ProbIniTreeInfection) * (1 - ADBInf);
								double noEABnoADB = (1 - ProbIniTreeInfection) * (1 - ADBInf);
								double noEABwithADB = (1 - ProbIniTreeInfection) * ADBInf;
								SetCrop(jcount, icount, 0, noEABnoADB); //All healthy
								SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
								SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
								SetCrop(jcount, icount, 3, noEABwithADB); //ADB: symp
								SetCrop(jcount, icount, 4, EABNotonADB); //EAB present first year + No ADB
								SetCrop(jcount, icount, 5, 0.0); //EAB present first year +ADB: Exposed not syp
								SetCrop(jcount, icount, 6, 0.0); //EAB present first year +ADB: Infected 
								SetCrop(jcount, icount, 7, EABonADB); //EAB present first year +ADB: symp
								SetCrop(jcount, icount, 8, 0.0); //EAB present longer + No ADB
								SetCrop(jcount, icount, 9, 0.0); //EAB present longer +ADB: Exposed not syp
								SetCrop(jcount, icount, 10, 0.0); //EAB present longer +ADB: Infected 
								SetCrop(jcount, icount, 11, 0.0); //EAB present longer +ADB: symp
							}
							//if (Sick==155000) // we plop in a seed of EAB of 1% of the trees
					
						}
						
							
					}
					
				}
				InF.close();
				InFF.close();
				InFFF.close();
				InF4.close();
				//InC.close();
				

			}
			else
				throw std::logic_error("Error opening inpur files");

			Grid::SetTotalPlantasFilled(-1);
			int mySymp(0); //Set to zero is max syptomatic trees in any one cell in the plantaion is less than 0.05
			double myWorry(0); //start unworried
			//int myCHMA(0);
			
			for (int pcount=0; pcount<MaxPlantID; pcount++)
			{
				//int myCells[MaxNumPCells]; 
				std::vector<int> myCells;
				int countcell=0;
				for (int jcount=0; jcount<Nrows; jcount++)
				{
					for (int icount=0; icount<Ncols; icount++)
					{
						if (pcount+1==PlantaID[jcount][icount])
						{
							myCells.push_back(CoordToGridID(jcount, icount)); //gathering the co-ordinates of all cells in a plantation(grove)
							//myCHMA=CHMAID[jcount][icount];
							
						}
					}
				}
				if (myCells.size()<1)
					throw std::logic_error("Number of cells in a plantation is less that 1");

		
				
			}

		}
		
	} //end of function SetLand

	void ReSetProbDetection()
	{
		for (int jcount = 0; jcount < Nrows; jcount++)
		{
			for (int icount = 0; icount < Ncols; icount++)
			{
				ProbDetection[jcount][icount] = 0;

			}
		}


	}
	void SetProbDetection()
	{
		int countcell = 0;
		for (int irow = 0; irow < Nrows; irow++)
		{
			for (int jcol = 0; jcol < Ncols; jcol++)
			{
				double myTreeDensity = Grid::GetTreeDensity(irow, jcol);
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
					double propEInfGTYr = tempTreeStates[8] + tempTreeStates[9] + tempTreeStates[10] + tempTreeStates[11]; // GT a year of EAB
					double NumTreesHealthy = Grid::GetTreeDensity(irow, jcol) * propENoEAB;  // NB Tree Density is actually number !
					double NumTreesInfested = Grid::GetTreeDensity(irow, jcol) * (1 - propENoEAB);
				
					double myNumLarvae(0);
			
					double c_parameter(0.026); // Girdelled 0.5096; Traps 0.052; //visual = 0.026 VASTHI LOOK HERE
					double NumTreeInspect(10);
					Grid::GetOWLav(irow, jcol, myNumLarvae);
					double phloemPerTree = Grid::GetPhloemPerTree();
					double myDetect = 0;
				
					///////////////////////////////////////////////////////////////////////
					double P_LarvToAdult(0.0);
					Grid::GetPLarv(irow, jcol, P_LarvToAdult);
					///////////////////////////////////////////////////////////////////////

						double AveLavPerTree = myNumLarvae * phloemPerTree * NumTreesInfested / (myTreeDensity* P_LarvToAdult);
						myDetect = 1 - exp(-c_parameter * NumTreeInspect * AveLavPerTree);
				
					if (myDetect != 0)
						int junk = 1;
					ProbDetection[irow][jcol] = ProbDetection[irow][jcol] + myDetect; //Note we add this up so very possible to be > 1
				}
				//else
				//{
					//ProbDetection[irow][jcol] = 0;
				//}

			}
		}

	}

	void SetWeightedProbDetection(double myWeight) // See document 
	{
		int countcell = 0;
		for (int irow = 0; irow < Nrows; irow++)
		{
			for (int jcol = 0; jcol < Ncols; jcol++)
			{
				double myTreeDensity = Grid::GetTreeDensity(irow, jcol);
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
					myDetect = 1 - exp(-c_parameter * NumTreeInspect * AveLavPerTree);

					if (myDetect != 0)
						int junk = 1;
					ProbDetection[irow][jcol] = ProbDetection[irow][jcol] + myWeight*myDetect; //Note we add this up so very possible to be > 1
				}
				/*else
				{
					ProbDetection[irow][jcol] = 0;
				}*/

			}
		}

	}

	double GetProbDetect(int irow, int jcol)
	{
		if ((irow > -1) && (irow < Grid::GetNumRows()) && (jcol > -1) && (jcol < Grid::GetNumCols()))
			return ProbDetection[irow][jcol];
		else
			throw std::logic_error("index out of bounds in GetProbDetect");

	}

	void SetLandBasic(int& iflag)
	{

		int UBlockID = 1; //unique Block ID	

		for (int icount = 0; icount < Nrows; icount++)
		{
			for (int jcount = 0; jcount < Ncols; jcount++)
			{
				SetTreeDensity(icount, jcount, 0.0); // This is just cleaning it before we start was SetCHMA(icount,jcount,0.0); 
				SetInitialTreeDensity(icount, jcount, 0.0); // This is just cleaning it before we start was SetCHMA(icount,jcount,0.0); 
				SetPlanta(icount, jcount, -9);
				//SetBlock(icount,jcount,-9); //removed all block stuff for now don't think we need it
				for (int kcount = 0; kcount < NStatus; kcount++)
					SetCrop(icount, jcount, kcount, 0);  //set all status to zero

			}
		}

		char path_buffer[_MAX_PATH];
		char Npath_buffer[_MAX_PATH];
		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);

		strcat_s(Npath_buffer, "\\TreeDensity\\Owners_GB.txt");

		std::ifstream    InF(Npath_buffer);

		strcpy_s(Npath_buffer, path_buffer);
		strcat_s(Npath_buffer, "\\TreeDensity\\ADBInfection_GB.txt");
		std::ifstream    InFF(Npath_buffer);

		strcpy_s(Npath_buffer, path_buffer);
		strcat_s(Npath_buffer, "\\TreeDensity\\treeDensity_GB.txt");
		std::ifstream	 InFFF(Npath_buffer);

		strcpy_s(Npath_buffer, path_buffer);
		strcat_s(Npath_buffer, "\\TreeDensity\\InclusionProbGB.txt");

		std::ifstream	 InF4(Npath_buffer);

		int TempOwner;
		double TempTreeDensity, ADBInf, infectionStatus;
		int Sick(0); //Flag to see of sick tree set/
		int iyr = 2;
		int MaxPlantID = -1; // This keeps track of the maximum ID number for the plantations for use later 
		
		///if (InF && InFF && InFFF)
		if (InF && InFF && InFFF && InF4)
		{
			double myCountProb = 0;
			int infCelliCol = 0;
			int infCelljRow = 0;
			for (int jcount = 0; jcount < Nrows; jcount++)
			{
				for (int icount = 0; icount < Ncols; icount++)
				{
					InF >> TempOwner;
					PlantaID[jcount][icount] = TempOwner;
					if (TempOwner > MaxPlantID)
						MaxPlantID = TempOwner;
					InFF >> ADBInf; // This is assumed to be symtomatic we make no attempt to predict other states
					InFFF >> TempTreeDensity;
					TempTreeDensity = TempTreeDensity*TreesPerHa;
					G_TreeDensity[jcount][icount] = TempTreeDensity;
					IniG_TreeDensity[jcount][icount] = TempTreeDensity;
					if (TempTreeDensity < 1.0)
						TempTreeDensity = 1;
					double ProbIniTreeInfection = 1 / TempTreeDensity; // This is now set so we have 1 tree infected per cell.

					double IncPr(0);
					InF4 >> IncPr; 
					SetInclustionProb(jcount, icount, IncPr/100.0);

					//Here we set the first cell that has a crop to 10% EAB landed 
					if (TempTreeDensity > 0)
					{
						SetCrop(jcount, icount, 0, 1 - ADBInf); //All healthy
						SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
						SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
						SetCrop(jcount, icount, 3, ADBInf); //ADB: symp
						SetCrop(jcount, icount, 4, 0.0); //EAB present not detected + No ADB
						SetCrop(jcount, icount, 5, 0.0); //EAB present not detected +ADB: Exposed not syp
						SetCrop(jcount, icount, 6, 0.0); //EAB present not detected +ADB: Infected 
						SetCrop(jcount, icount, 7, 0.0); //EAB present not detected +ADB: symp
						SetCrop(jcount, icount, 8, 0.0); //EAB detected + No ADB
						SetCrop(jcount, icount, 9, 0.0); //EAB detected +ADB: Exposed not syp
						SetCrop(jcount, icount, 10, 0.0); //EAB detected +ADB: Infected 
						SetCrop(jcount, icount, 11, 0.0); //EAB detected +ADB: symp

						if ((jcount == infCelljRow) && (icount == infCelliCol))
						{
							double EABonADB = (ProbIniTreeInfection)*ADBInf;
							double EABNotonADB = (ProbIniTreeInfection) * (1 - ADBInf);
							double noEABnoADB = (1 - ProbIniTreeInfection) * (1 - ADBInf);
							double noEABwithADB = (1 - ProbIniTreeInfection) * ADBInf;
							SetCrop(jcount, icount, 0, noEABnoADB); //All healthy
							SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
							SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
							SetCrop(jcount, icount, 3, noEABwithADB); //ADB: symp
							SetCrop(jcount, icount, 4, EABNotonADB); //EAB present first year + No ADB
							SetCrop(jcount, icount, 5, 0.0); //EAB present first year +ADB: Exposed not syp
							SetCrop(jcount, icount, 6, 0.0); //EAB present first year +ADB: Infected 
							SetCrop(jcount, icount, 7, EABonADB); //EAB present first year +ADB: symp
							SetCrop(jcount, icount, 8, 0.0); //EAB present longer + No ADB
							SetCrop(jcount, icount, 9, 0.0); //EAB present longer +ADB: Exposed not syp
							SetCrop(jcount, icount, 10, 0.0); //EAB present longer +ADB: Infected 
							SetCrop(jcount, icount, 11, 0.0); //EAB present longer +ADB: symp
						}

					}


				}

			}
			InF.close();
			InFF.close();
			InFFF.close();
			InF4.close();
			//InC.close();


		}
		else
			throw std::logic_error("Error opening inpur files");

		double P_LarvToAdultStart = Grid::GetPropOverWinterLarvaeStart();

		for (int icount = 0; icount < Grid::GetNumRows(); icount++)
		{
			for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
			{
				if ((Grid::GetCrop(icount, jcount, 4) + Grid::GetCrop(icount, jcount, 7)) > 0) ////Status of trees. Status 4 is 'removed which includes never planted there. If this is the case there will be no bug
				{


					Grid::SetNumG(icount, jcount, 10.0); //Sets number healthy bugs
					Grid::SetOWLav(icount, jcount, 0);
					Grid::SetTempBug(icount, jcount, 0);
					GetPropOverWinterLarvaeStart();
					Grid::SetInvasion(icount, jcount, -1);


				}
				else
				{
					Grid::SetNumG(icount, jcount, 0); //Sets number healthy bugs
					Grid::SetOWLav(icount, jcount, 0);
					Grid::SetTempBug(icount, jcount, 0);
					Grid::SetPLarv(icount, jcount, P_LarvToAdultStart);
					Grid::SetInvasion(icount, jcount, -1);
				}


			}
		}

		Grid::SetTotalPlantasFilled(-1);
		int mySymp(0); //Set to zero is max syptomatic trees in any one cell in the plantaion is less than 0.05
		double myWorry(0); //start unworried
		//int myCHMA(0);




	} //end of function SetLandBasic

	void ReSetLandBasic( int infCelljRow, int infCelliCol, int& iflag) 
	{

		int UBlockID = 1; //unique Block ID	
		char path_buffer[_MAX_PATH];
		char Npath_buffer[_MAX_PATH];
		_getcwd(path_buffer, _MAX_PATH);
		
		strcpy_s(Npath_buffer, path_buffer);
		//strcat_s(Npath_buffer, "\\TreeDensity\\UKGridADB.txt"); 
		//strcat_s(Npath_buffer, "\\TreeDensity\\ADBInfection_SE.txt"); 
		strcat_s(Npath_buffer, "\\TreeDensity\\ADBInfection_GB.txt");
		std::ifstream    InFF(Npath_buffer);
		double ADBInf(0);
		
		
		for (int jcount = 0; jcount < Nrows; jcount++)
		{
			for (int icount = 0; icount < Ncols; icount++)
			{
				double TempTreeDensity = Grid::GetInitialTreeDensity(jcount, icount);
				if (TempTreeDensity < 1.0)
					TempTreeDensity = 1;
				double ProbIniTreeInfection = 1 / TempTreeDensity; // This is now set so we have 1 tree infected per cell.
				InFF >> ADBInf; // This is assumed to be symtomatic we make no attempt to predict other states
				ADBInf = 0;
					//Here we set the first cell that has a crop to 10% EAB landed 
					if (TempTreeDensity > 0)
					{
						Grid::SetTreeDensity(jcount, icount, Grid::GetInitialTreeDensity(jcount, icount)); // NATHAN THIS IS THE FIX
						SetCrop(jcount, icount, 0, 1 - ADBInf); //All healthy
						SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
						SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
						SetCrop(jcount, icount, 3, ADBInf); //ADB: symp
						SetCrop(jcount, icount, 4, 0.0); //EAB present not detected + No ADB
						SetCrop(jcount, icount, 5, 0.0); //EAB present not detected +ADB: Exposed not syp
						SetCrop(jcount, icount, 6, 0.0); //EAB present not detected +ADB: Infected 
						SetCrop(jcount, icount, 7, 0.0); //EAB present not detected +ADB: symp
						SetCrop(jcount, icount, 8, 0.0); //EAB detected + No ADB
						SetCrop(jcount, icount, 9, 0.0); //EAB detected +ADB: Exposed not syp
						SetCrop(jcount, icount, 10, 0.0); //EAB detected +ADB: Infected 
						SetCrop(jcount, icount, 11, 0.0); //EAB detected +ADB: symp

						if ((jcount == infCelljRow) && (icount == infCelliCol))
						{
							double EABonADB = (ProbIniTreeInfection)*ADBInf;
							double EABNotonADB = (ProbIniTreeInfection) * (1 - ADBInf);
							double noEABnoADB = (1 - ProbIniTreeInfection) * (1 - ADBInf);
							double noEABwithADB = (1 - ProbIniTreeInfection) * ADBInf;
							SetCrop(jcount, icount, 0, noEABnoADB); //All healthy
							SetCrop(jcount, icount, 1, 0.0); // ADB: Exposed not syp
							SetCrop(jcount, icount, 2, 0.0); //ADB: Infected 
							SetCrop(jcount, icount, 3, noEABwithADB); //ADB: symp
							SetCrop(jcount, icount, 4, EABNotonADB); //EAB present first year + No ADB
							SetCrop(jcount, icount, 5, 0.0); //EAB present first year +ADB: Exposed not syp
							SetCrop(jcount, icount, 6, 0.0); //EAB present first year +ADB: Infected 
							SetCrop(jcount, icount, 7, EABonADB); //EAB present first year +ADB: symp
							SetCrop(jcount, icount, 8, 0.0); //EAB present longer + No ADB
							SetCrop(jcount, icount, 9, 0.0); //EAB present longer +ADB: Exposed not syp
							SetCrop(jcount, icount, 10, 0.0); //EAB present longer +ADB: Infected 
							SetCrop(jcount, icount, 11, 0.0); //EAB present longer +ADB: symp
						}

					}


				}

			}
			
		InFF.close();
		double P_LarvToAdultStart = GetPropOverWinterLarvaeStart();
		for (int icount = 0; icount < Grid::GetNumRows(); icount++)
		{
			for (int jcount = 0; jcount < Grid::GetNumCols(); jcount++)
			{
				if ((Grid::GetCrop(icount, jcount, 4) + Grid::GetCrop(icount, jcount, 7)) > 0) ////Status of trees. Status 4 is 'removed which includes never planted there. If this is the case there will be no bug
				{


					Grid::SetNumG(icount, jcount, 10.0); //Sets number healthy bugs
					Grid::SetOWLav(icount, jcount, 0);
					Grid::SetTempBug(icount, jcount, 0);
					Grid::SetPLarv(icount, jcount, P_LarvToAdultStart);
					Grid::SetInvasion(icount, jcount, -1);


				}
				else
				{
					Grid::SetNumG(icount, jcount, 0); //Sets number healthy bugs
					Grid::SetOWLav(icount, jcount, 0);
					Grid::SetTempBug(icount, jcount, 0);
					Grid::SetPLarv(icount, jcount, P_LarvToAdultStart);
					Grid::SetInvasion(icount, jcount, -1);
				}


			}
		}


	} //end of function SetLandBasic
}


////////////////////////////END OF NAMESPACE GRID///////////////////////////////////////////

namespace FlightP // keeps a vector of the proportions of Bugs of each type that fly each distance
	
{
	std::vector<double> BugFlight; //
	
	void ClearData()
	{
		
		BugFlight.clear(); //
		

	};


	void MakeProportions(std::vector<double>& Type, double lambda)
	{
			//Type 1 (Type[0]) centre cell for other types see 
			int idist=0;
			int jdist=0;
			
			Type.push_back(Integrate(idist, jdist, lambda)); //position 0
			Type.push_back(Integrate(0, 1,  lambda)); //position 1
			Type.push_back(Integrate(1, 1, lambda)); //position 2
			//Test to see how much of the population is in this inner square
			double popT=Type[0]+4*Type[1]+4*Type[2];
			//if (popT<0.99) //keep calculating
			int Mcount=3;
			while (popT<0.95) 
			{
				for (int icount=0; icount<Mcount; icount++)
				{
					double Ival=Integrate(icount, Mcount-1, lambda);
					Type.push_back(Ival);
					if ((icount>0)&&(icount<Mcount-1))
						popT=popT+8*Ival;
					else
						popT=popT+4*Ival;
				}
				Mcount++;
				if (Mcount>900)
					int junk=0;
								
			}
			//put remained in the centre
			
			Type[0]=Type[0]+1.0-popT;
			if (Type[0]<0) Type[0]=0;
	

	}

	void Initialise()
	{ 
		//mercader 2009 refit o.o12 so 12 if per km. // mercader 2016 0.0017, so 1.7 // exqamples i used 4 - guess!
		//Scrader - EU experts Guessed 50% chance within 1.5km - i think this is after detection?? otherwise kernel very wide
		// for Cauchy using mercader 2009 25.5m so lambda = 0.0255
		double lambda = 0.0414; // 0.0255; // 0.00035; // For weekly time step. The value of 0.46 is based on Schrader et al. 2020 50% traveled 1.5 km in a year
			MakeProportions(BugFlight, lambda);
		


	}

	void GetBugFlight(std::vector<double>& MyDist)
	{
		int Len=BugFlight.size();
		for (int icount=0; icount<Len; icount++)
		{
			MyDist.push_back(BugFlight[icount]);
		}
	}

	


		
} //end of namespace

double Integrate(int idist, int jdist, double lambda)
{
	//This function integrates the exponential distributions over each cell to see how many bugs we have in a cell.

	//idist= irow-i0 where irow is place we are estimating population and i0 is where Buge started
	//jdist=jrow-j0 ditto
	double MygridLen = Grid::GetGridLen();
	
	Grid::SetIntVars(jdist, MygridLen, lambda);

	double xlb[2], xub[2];
	xlb[0]=(idist -0.5)* MygridLen;
	xlb[1]=(jdist -0.5)* MygridLen;
	xub[0]=(idist +0.5)* MygridLen;
	xub[1]=(jdist +0.5)* MygridLen;

	

	
	double ABSACC=0.0;
	double RELACC=0.0000001;
	//double ANS(0);
	int NDIM(2);
	int NumFunc(1);
	int MaxCals(18000);
	int MinCals(10);
	const int LENWRK=50000;//(NDIM+NumFunc+2)*(10+MaxCals);
	double WRKSTR[LENWRK];
	int IFAIL=1;
	double ANS[1];
	double ABSET[1];
		
	
	//D01DAF(xlb, xub, ylbFunc, yubFunc,fFunc,ABSACC, ANS, NPTS, IFAIL);
	D01EAF(NDIM, xlb, xub, MinCals, MaxCals, NumFunc, fFunc2, ABSACC, RELACC, LENWRK, WRKSTR, ANS, ABSET, IFAIL);

	if (ANS[0]<0)
		throw std::logic_error("Error in integration");
	if (IFAIL!=0)
		int junk=0;
	
	//ANS[0]=ANS[0]*ga;
	return ANS[0];

}


extern "C" void __stdcall fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] )
//void  fFunc2(const int& NDIM, const double x[], const int& NumFunc, double funcY[] )
{
	double jdist, gridLen, lambda;
	Grid::GetIntVars(jdist, gridLen, lambda);
	int IFAIL=1;
	double r=pow(x[0],2.0)+pow(x[1],2.0);
	r=pow(r, 0.5);
	double func(0);
	double OOtwoPi=1.0/(2.0*3.141572);
	if (r>0)
	{
		//func=OOtwoPi*lambda*exp(-lambda*r); // exponential pdf divide by circumference // updated but check
		// Cauchy some issue rescaling to km // need to understand r 
		func = lambda * lambda / ((3.14159265359 * 3.14159265359) * (x[0] * x[0] + lambda * lambda) * (x[1] * x[1] + lambda * lambda)); // Cauchy pdf divide by circumference
	}
	else
	{
		func=0;
	}

	funcY[0]=func;


}

void GridIDToCoord(int Floc, int& irow, int& icol)
{
	irow = Floc % Grid::GetNumRows();  // % is remainder after dividing
	icol=(floor(Floc/Grid::GetNumRows())); 
}
int CoordToGridID(int jrow, int icol)
{
	return icol*Grid::GetNumRows()+jrow;
}



