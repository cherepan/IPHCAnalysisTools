#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];
$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";
if($UserID eq "vcherepa"){
    $UserIDCern="cherepan";
}
$CMSSWRel = "8_0_25";
$ARCH='slc6_amd64_gcc530';
printf("\nYour user ID is:   $UserID");
if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\nThe installation follows the recomendation given in the description");
    printf("\nof LLRHiggsTauTau for the latest CMSSW release: https://github.com/LLRCMS/LLRHiggsTauTau");
    printf("\n\nThis code requires one input option. The syntax is: ./todo.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --TauNtuple <TauNtupleDir>               Setups up TauNtuple and gives instructions for submission");
    printf("\n                                                   <TauNtupleDir> location of CMSSW for TauNtuple.");
    printf("\n                                                      Optional Commmads: ");
    printf("\n                                                     --CMSSWRel <CMSSW release #> The CMSSW release you want to use Default: $CMSSWRel");
    printf("\n                                                     --PileupVersion <Version # ie V08-03-17>");
    printf("\n                                                     --BTag <YES/NO> Default: Yes (To be implemented)");
    printf("\n                                                     --SVfit branch (To be implemented). Default is svFit_2015Apr03");
    printf("\n                                                     --SVfit Option to turn on SVfit  (To be implemented) \n");
    exit(0);  
}
 my $dir = getcwd;
$time= strftime("%h_%d_%Y",localtime);

if( $ARGV[0] eq "--TauNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --TauNtuple <TauNtupleDir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }
    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir \n");

    system(sprintf("rm Install_TauNtuple_$time"));

    system(sprintf("echo \"cernlib-use root\" >> Install_TauNtuple_$time"));

    system(sprintf("echo \"export SCRAM_ARCH=\\\"$ARCH\\\"\" >> Install_TauNtuple_$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_TauNtuple_$time"));


   system(sprintf("echo \"cd $basedir\" >>  Install_TauNtuple_$time")); 
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >>  Install_TauNtuple_$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_TauNtuple_$time")); 
    system(sprintf("echo \"cmsenv\" >> Install_TauNtuple_$time")); 


   system(sprintf("echo \"git cms-merge-topic cms-met:METRecipe_8020\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git cms-merge-topic ikrav:egm_id_80X_v2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git cms-merge-topic gpetruc:badMuonFilters_80X_v2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/LLRCMS/LLRHiggsTauTau\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd LLRHiggsTauTau; git checkout master\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/latinos/UserCode-sixie-Muon-MuonAnalysisTools Muon/MuonAnalysisTools\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd Muon/MuonAnalysisTools ; git checkout master -- interface/MuonEffectiveArea.h\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/cms-analysis/EgammaAnalysis-ElectronTools EGamma/EGammaAnalysisTools\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd EGamma/EGammaAnalysisTools; git checkout c0db796 -- interface/ElectronEffectiveArea.h\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone -n https://github.com/VBF-HZZ/UFHZZAnalysisRun2\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd UFHZZAnalysisRun2 ; git checkout master FSRPhotons\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd -\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/veelken/SVfit_standalone TauAnalysis/SVfitStandalone\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd TauAnalysis/SVfitStandalone\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git checkout svFit_2015Apr03\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/src\" >> Install_TauNtuple_$time")); 

   system(sprintf("echo \"scram b -j 4\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/external/$SCRAM_ARCH\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd data/RecoEgamma/ElectronIdentification/data\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"git checkout egm_id_80X_v1\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"cd $CMSSW_BASE/src\" >> Install_TauNtuple_$time")); 
   system(sprintf("echo \"scram b -j 4\" >> Install_TauNtuple_$time")); 
    # print Instructions
    printf("\n\nInstructions:");
    printf("\nsource  Install_TauNtuple_$time to complete installation, compilation might take some time...  \n\n");
    printf("\nTo run test job do  'cmsRun analyzer.py'  in  $CMSSW_BASE/src/LLRHiggsTauTau/NtupleProducer/test \n\n");

}
