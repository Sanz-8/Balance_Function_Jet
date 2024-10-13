#include "TString.h"
#include <iostream> 
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include <fstream>

void File_Write(bool isTestRun = false, string inputfile = "Step_1.txt", string outfile = "Step1", string save_outrootfile_folder = "/eos/user/r/rpradhan/Research/With_NirbhayBhai", int nCpus = 4, string jobflavour= "workday")
{
  std::string voms = "no";
  std::cout<<"did you do voms (yes or no)? if not do that first"<<std::endl;
  cin>>voms;
  if(voms == "no")
    {
      gSystem->Exec("voms-proxy-init -rfc -voms cms --out voms_proxy.txt --hours 200");
    }
  std::cout<<std::endl;

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

	  if(isTestRun && nlines >= 5) 
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

      gSystem->Exec(Form("mkdir -p %s/%s", save_outrootfile_folder.c_str(), outfile.c_str())); //create the folder to save the output root file

      gSystem->Exec("mkdir -p condor_out"); // create the folder to save .log, .err, .out file                                      
      gSystem->Exec(Form("mkdir -p condor_out/%s", outfile.c_str())); // create the folder to save .log, .err, .out file             
      gSystem->Exec(Form("rm -rf condor_out/%s/*", outfile.c_str())); // delete existing things 
      gSystem->Exec("mkdir -p Compiler_Files"); //create a folder to contain Compiler_Files                                      
      gSystem->Exec(Form("mkdir -p Compiler_Files/%s", outfile.c_str())); //create a folder to contain Submission_Files             
      gSystem->Exec(Form("rm -rf Compiler_Files/%s/*", outfile.c_str())); // delete existing things  
      gSystem->Exec(Form("g++ -o Compiler_Files/%s/compile_%s.out Tree_Analyzer_NewNew_Condor_CMSpythia.C `root-config --cflags --glibs`", outfile.c_str(),outfile.c_str()));
      
      gSystem->Exec("mkdir -p Executable_Files"); //create a folder to contain Submission_Files
      gSystem->Exec(Form("mkdir -p Executable_Files/%s", outfile.c_str())); //create a folder to contain Submission_Files        
      gSystem->Exec(Form("rm -rf Executable_Files/%s/*", outfile.c_str())); // delete existing things 

      ofstream myfile_exe;
      myfile_exe.open(Form("Executable_Files/%s/submit_%s.sh", outfile.c_str(), outfile.c_str()));
      myfile_exe<<"#!/bin/bash"<<std::endl;
      myfile_exe<<std::endl;
      myfile_exe<<"cd /afs/cern.ch/user/s/slaishra/project/Injet_BalanceFUnction_Pythia/CMS_Data_pp_13TeV_parker"<<std::endl;
      myfile_exe<<"echo \"Submit skim jobs at \""<<std::endl;
      myfile_exe<<"echo PWD: $PWD"<<std::endl;
      myfile_exe<<""<<std::endl;
      myfile_exe<<"./Compiler_Files/"<<outfile.c_str()<<"/"<<"compile_"<<outfile.c_str()<<".out"<<" "<<"\"$1\" \"$2\" \"$3\" \"$4\""<<std::endl;
      myfile_exe<<endl;
      myfile_exe.close();
      
      gSystem->Exec(Form("chmod u+x Executable_Files/%s/submit_%s.sh", outfile.c_str(), outfile.c_str()));
      gSystem->Exec("mkdir -p Submission_Files"); //create a folder to contain Submission_Files                                      
      gSystem->Exec(Form("mkdir -p Submission_Files/%s", outfile.c_str())); //create a folder to contain Submission_Files  
      gSystem->Exec(Form("rm -rf Submission_Files/%s/*",outfile.c_str())); // delete existing things 
      
      ofstream myfile;
      myfile.open(Form("Submission_Files/%s/Submission_%s.sub", outfile.c_str(), outfile.c_str()));
      myfile<<"universe     = vanilla"<<endl;
      myfile<<"getenv       = True"<<endl;
      myfile<<"executable   = Executable_Files/"<<outfile.c_str()<<"/"<<"submit_"<<outfile.c_str()<<".sh"<<endl;
      myfile<<"+JobFlavour  = \""<<jobflavour.c_str()<<"\""<<endl;
      myfile<<"requirements = (OpSysAndVer =?= \"CentOS7\")"<<endl;
      myfile<<"RequestCpus  = "<<nCpus<<endl;
      myfile<<"transfer_input_files  = voms_proxy.txt"<<endl;
      myfile<<"environment = \"X509_USER_PROXY=voms_proxy.txt\""<<endl;
      myfile<<endl;
      myfile<<"log          = condor_out/"<<outfile.c_str()<<"/"<<outfile.c_str()<<"_$(process)"<<".log"<<endl;
      myfile<<"output          = condor_out/"<<outfile.c_str()<<"/"<<outfile.c_str()<<"_$(process)"<<".out"<<endl;
      myfile<<"error          = condor_out/"<<outfile.c_str()<<"/"<<outfile.c_str()<<"_$(process)"<<".err"<<endl;
      myfile<<"arguments    = Divided_txt_files/"<<outfile.c_str()<<"/Outfile_"<<outfile.c_str()<<"_$(process)"<<".txt"<<"  "<<"$(process)"<<"  "<<save_outrootfile_folder.c_str()<<"/"<<outfile.c_str()<<"  "<<outfile.c_str()<<endl;
      myfile<<"queue "<<numJobs<<endl;
      myfile<<endl;
      myfile.close();
    }
}
