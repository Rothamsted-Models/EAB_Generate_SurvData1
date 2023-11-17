#ifndef __Model_H__
#define __Model_H__

//Constants for optimisation
const int ConstDefNumLocs(5); // the number of locations that we might sample at
void InitialiseModel(); //  proportion that believe control works
void SetLand(int readin); // 0 to read from file and 1 to create, proportion that believe control works
void DeleteModel();
void ModelMonth(int imth, double randInf, double kill_rate); //if randInf=0 then no random infection event if 1 then a randomly chosed cell has a single infected bug land with probability randIf
void Generation(int irow, int jcol, int iyr, double kill_rate);
void Tree_Generation(int irow, int jcol, int iyr); //progesses tree health status in a time step
void Flight();
void BugsJump();
void Reflect(int& nrow, int& ncol);
void Distribute2(int irow, int jcol, double num);
void flyAgain(int& nrow, int& ncol);
void AdjustEABStatusPostFlight(int irow, int jcol, int FellingType, double propFelled);


void WriteCrops(int imth); //writes crops for a particular year
void WriteLarvT(int imth); //writes total Larvae
void WriteBugsT(int imth); //writes total bugs 
void DetectionFullRunW_Target(double RandIf, double kill_rate,  int ConstMaxInvYr);
void WriteAveLarvae(int imth, int Repeat); //writes LArvae per tree for detection for a particular year
void WritePropDead(int imth, int Repeat);
void WriteProbDetectForR();
void SaveDataInStruc(double SaveStructure[], double SaveStructureT[], int tcount);
void WriteDetectionInforOpt_Track(int occassion, int infRow, int infCol, double myWeight, int MaxYrs, double SaveStructure[], double SaveStructureT[]);
#endif