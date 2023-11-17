//This code was developed by Alice Milne,  Vasthi Alonso Chavez and Nathan Brown Rothamsted Research, Harpenden, UK
//Please contact alice.milne@rothamsted.ac.uk with any questions
//The code was developed as part of the NERC funded Smarties project RP10560-10


//*********************************************Purpose of code*************************************
//The code simulates the invasion and spread of EAB across GB. Details are given in Brown et al., 2023
//This code uses an input file "SampleEntry.txt" to define the initial location for about 10000 runs. 
// This file was generated with the associated Matlab code that defines risk based on firewood imports and use.
//The outputs are stored in DataForOpt.txt. This output file is then used in a seperate optimisation step to define optimised sampling.
//The reason for decoupling the main code from optimisation was speed and flexibility.
//We have attenpted to clean up the code so as not to confuse the reader, but there are likely to be relics from other versions and other applications of this code.
//We hope you have as much fun with it as we did.



// EAB_V1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <stdexcept>
#include  <fstream>
#include <vector>
#include <direct.h>
#include "EABModel.h"
#include "EABGrid.h"

int main()
{

	//Initialise the landscape
	try
	{
		int iniMode = 1; //0 = uniform trees no landscape structure; 1= read from files; 2= uniform but has plantations and blocks
		double RandIf = 0.0; // The probability that a single bug arrives at a randomly chosen host cell
		
		//Control parameters
		double kill_rate = 0.94; //Not used in this model but there for future
		const int ConstMaxInvYr(8); // The number of years we look to "find" EAB for the optimisation
		
		//Set up some write files
		char path_buffer[_MAX_PATH];
		char Npath_buffer[_MAX_PATH];
		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);
		

		InitialiseModel( ); // initialise the model
		
		//WriteCrops(0); // This writes out the various infection states so will need revising when we know what we want to write out. Commented out for now
		
		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);
		strcat_s(Npath_buffer, "\\OutFiles\\MyErrors.txt");
		std::ofstream of(Npath_buffer);
		of.close();

		_getcwd(path_buffer, _MAX_PATH);
		strcpy_s(Npath_buffer, path_buffer);
		char ECBStr[50];
		strcpy_s(ECBStr, "\\OutFiles\\DataForOpt.txt");
		strcat_s(Npath_buffer, ECBStr);
		std::ofstream OpOutF(Npath_buffer);
		OpOutF.close();

		Grid::ReSetProbDetection();
		int iflag(0);
		Grid::SetLandBasic(iflag); // Very basic setup just for this optimisation 
		DetectionFullRunW_Target(RandIf, kill_rate,  ConstMaxInvYr); // This uses preselected start points
		
		std::cout << "Job complete";


	}
	catch (std::logic_error& E)
	{
		char Mess[500];
		strcpy_s(Mess, E.what());
		std::cout << Mess;


	}
    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
