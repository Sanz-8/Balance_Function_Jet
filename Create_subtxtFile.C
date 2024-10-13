#include "TString.h"
#include <iostream> 
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include <fstream>

void Create_subtxtFile(bool isTestRun = false, string inputfile = "Step_1.txt", string outfile = "Step1")
{
  int nlines = 0; 
  ifstream fpr(inputfile.c_str(), ios::in);

  string line;

  std::vector<string> file_name_vector;
  
  if(fpr.is_open()) 
    {
      while(fpr.peek() != EOF)
	{
	  getline(fpr, line);
	  nlines++;

	  if(isTestRun && nlines > 400) 
	    {
	      std::cout<<"only "<<nlines-1<<" root files are stored for test run"<<std::endl;
	      break;
	    }
 
	  file_name_vector.push_back(line);
	}
    }

  fpr.close();

  std::cout<<endl;
  std::cout<<"total file size is: "<<file_name_vector.size()<<std::endl;
  std::cout<<endl;
  // try to give an input value for njobs such that, it would be ((njobs*ratio)-file_name_vector.size()) < ratio and files will approximately divided equally.

  int njobs;

  if(isTestRun)
    {
      std::cout<<"~~~~~~~~~~~~:Enter the njobs for test run:~~~~~~~~~~~~~~~~~"<<std::endl;
      std::cout<<"Enter the number of jobs you want: "<<std::endl; 
    }
  else
    {
      std::cout<<"~~~~~~~~~~~~:Enter the njobs for final run:~~~~~~~~~~~~~~~~~"<<std::endl;
      std::cout<<"Enter the number of jobs you want: "<<std::endl;
    }
  
  cin>> njobs;

  int ratio = file_name_vector.size()/njobs; 

  ofstream fout;
  
  if(ratio < 1)
    {
      std::cout<<"Number of jobs greater than number of files"<<std::endl;
      std::cout<<"Number of jobs should be less than number of files, please reduce the number of jobs"<<std::endl;
      exit(0);
    }

  if(ratio > file_name_vector.size())
    {
      std::cout<<"it will never happen"<<std::endl;
      exit(0);
    }

  if(ratio != file_name_vector.size() || ratio == file_name_vector.size())
    {
      int numJobs = njobs;
      
      gSystem->Exec("mkdir -p Divided_txt_files"); //create the folder to save the divided .txt files
      gSystem->Exec(Form("mkdir -p Divided_txt_files/%s", outfile.c_str())); //create the folder to save the divided .txt files
      gSystem->Exec(Form("rm -rf Divided_txt_files/%s/*", outfile.c_str())); // delete previously created .txt files
            
      for(int i = 0; i < njobs; i++)
	{
	  int start = i*ratio;
	  int end = (i + 1) * ratio;
	  
	  if(i == njobs - 1)
	    {
	      start = i*ratio;
	      end = file_name_vector.size();
	    }
	  
	  fout.open(Form("Divided_txt_files/%s/Outfile_%s_%d.txt", outfile.c_str(), outfile.c_str(), i));
	  
	  for(int k = start; k<end; k++)
	    {
	      fout <<file_name_vector[k].c_str()<<std::endl;
	    }
	  fout.close();
	}
    }
}
