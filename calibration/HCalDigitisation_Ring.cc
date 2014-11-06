/**
 *  @file   PandoraAnalysis/calibration/HCalDigitisation_Ring.cc
 * 
 *  @brief  Mean of distribution of direction corrections applied to sim calorimeter hits in HCal Barrel, EndCap and Other 
 *          for HCal ring digitisation
 * 
 *  $Log: $
 */
#include "TApplication.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1F.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TTree.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdio.h>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

/**
 *  @brief  Parameters class
 */
class HCalRingDigitisation 
{
public:
    /**
     *  @brief  Constructor
    */
    HCalRingDigitisation();

    /**
     *  @brief  Destructor
     */
    ~HCalRingDigitisation();

    /**
     *  @brief  Create histograms
    */
    void MakeHistograms();

// Non trivial setting on initialisation
    float           m_trueEnergy;                           ///< True energy (opposed to kinetic) of particle being simulated
    std::string     m_inputMuonRootFiles;                   ///< Input root files
    std::string     m_outputPath;                           ///< Output path to send results
    float           m_PeakHCalBarrelDirectionCorrectedADC;  ///< Peak position in m_hHCalBarrelDirectionCorrectedADC
    float           m_PeakHCalEndCapDirectionCorrectedADC;  ///< Peak position in m_hHCalEndCapDirectionCorrectedADC
    float           m_PeakHCalOtherDirectionCorrectedADC;   ///< Peak position in m_hHCalOtherDirectionCorrectedADC

private:
    /**
     *  @brief  Histogram peak finder
    */
    int PeakFinder(const TH1F *const pTH1F);

    TH1F           *m_hHCalBarrelDirectionCorrectedADC;     ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal Barrel
    TH1F           *m_hHCalEndCapDirectionCorrectedADC;     ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal EndCap
    TH1F           *m_hHCalOtherDirectionCorrectedADC;      ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal Other
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 * 
 *  @return success
 */

bool ParseCommandLine(int argc, char *argv[], HCalRingDigitisation &hCalRingDigitisation);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        HCalRingDigitisation hCalRingDigitisation;

        if (!ParseCommandLine(argc, argv, hCalRingDigitisation))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        hCalRingDigitisation.MakeHistograms();

        std::string dataFileName(hCalRingDigitisation.m_outputPath + "Calibration.txt");
        ofstream    data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "_____________________________________________________________________________________" << std::endl;
        std::cout << "_____________________________________________________________________________________" << std::endl;

        data_file << "Digitisation of the HCal Ring                      : " << std::endl << std::endl;
        data_file << "For Muons with energy                              : " << hCalRingDigitisation.m_trueEnergy << " : " <<std::endl;
        data_file << "The MIP peaks for each part of the HCal were found : " << std::endl;
        data_file << "to occur at the following ADC values:              : " << std::endl;
        data_file << "HCal Barrel MIP Peak:                              : " << hCalRingDigitisation.m_PeakHCalBarrelDirectionCorrectedADC << " : " <<std::endl;
        data_file << "HCal EndCap MIP Peak:                              : " << hCalRingDigitisation.m_PeakHCalEndCapDirectionCorrectedADC << " : " <<std::endl;
        data_file << "HCal Ring MIP Peak:                                : " << hCalRingDigitisation.m_PeakHCalOtherDirectionCorrectedADC << " : " <<std::endl;

        std::cout << "Digitisation of the HCal Ring                      : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << hCalRingDigitisation.m_trueEnergy << " : " <<std::endl;
        std::cout << "The MIP peaks for each part of the HCal were found : " << std::endl;
        std::cout << "to occur at the following ADC values:              : " << std::endl;
        std::cout << "HCal Barrel MIP Peak:                              : " << hCalRingDigitisation.m_PeakHCalBarrelDirectionCorrectedADC << " : " <<std::endl;
        std::cout << "HCal EndCap MIP Peak:                              : " << hCalRingDigitisation.m_PeakHCalEndCapDirectionCorrectedADC << " : " <<std::endl;
        std::cout << "HCal Ring MIP Peak:                                : " << hCalRingDigitisation.m_PeakHCalOtherDirectionCorrectedADC << " : " <<std::endl;

        data_file.close();
        delete pTApplication;
    }

    catch (std::exception &exception)
    {
        std::cout << "Exception caught " << exception.what() << std::endl;
        delete pTApplication;
        return 1;
    }

    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

HCalRingDigitisation::HCalRingDigitisation() :
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_inputMuonRootFiles(""),
    m_outputPath(""),
    m_PeakHCalBarrelDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_PeakHCalEndCapDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_PeakHCalOtherDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_hHCalBarrelDirectionCorrectedADC(NULL),
    m_hHCalEndCapDirectionCorrectedADC(NULL),
    m_hHCalOtherDirectionCorrectedADC(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HCalRingDigitisation::~HCalRingDigitisation()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HCalRingDigitisation::MakeHistograms()
{
    m_hHCalBarrelDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCBarrel", "Distribution of Direction Corrected ADCs in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalBarrelDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalBarrelDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCEndCap", "Distribution of Direction Corrected ADCs in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalEndCapDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalEndCapDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectedADC = new TH1F("HCalDirectionCorrectedADCOther", "Distribution of Direction Corrected ADCs in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.005);
    m_hHCalOtherDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hHCalOtherDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

    unsigned int found_slash = m_inputMuonRootFiles.find_last_of("/");
    std::string path = m_inputMuonRootFiles.substr(0,found_slash);
    std::string file = m_inputMuonRootFiles.substr(found_slash+1);

    unsigned int found_star = file.find_last_of("*");
    std::string file_prefix = file.substr(0,found_star);
    std::string file_suffix = file.substr(found_star+1);

    TSystemDirectory dir(path.c_str(), path.c_str());
    TList *files = dir.GetListOfFiles();

    if (files)
    {
        TSystemFile *pTSystemFile;
        TString fname;
        TIter next(files);
        while ((pTSystemFile = (TSystemFile*)next()))
        {
            fname = pTSystemFile->GetName();

            if (!pTSystemFile->IsDirectory() && fname.EndsWith(file_suffix.c_str()) && fname.BeginsWith(file_prefix.c_str()))
            {
                TString filename(path + "/" + fname);
                TFile *i_file = new TFile(filename);
                TH1F *i_hHCalBarrelDirectionCorrectedADC = (TH1F*) i_file->Get("HCalDirectionCorrectedADCBarrel");
                TH1F *i_hHCalEndCapDirectionCorrectedADC = (TH1F*) i_file->Get("HCalDirectionCorrectedADCEndCap");
                TH1F *i_hHCalOtherDirectionCorrectionADC = (TH1F*) i_file->Get("HCalDirectionCorrectedADCOther");

                if (i_hHCalBarrelDirectionCorrectedADC!=NULL)
                    m_hHCalBarrelDirectionCorrectedADC->Add(i_hHCalBarrelDirectionCorrectedADC,1.0);

                if (i_hHCalEndCapDirectionCorrectedADC!=NULL)
                    m_hHCalEndCapDirectionCorrectedADC->Add(i_hHCalEndCapDirectionCorrectedADC,1.0);

                if (i_hHCalOtherDirectionCorrectionADC!=NULL)
                    m_hHCalOtherDirectionCorrectedADC->Add(i_hHCalOtherDirectionCorrectionADC,1.0);

                delete i_hHCalBarrelDirectionCorrectedADC;
                delete i_hHCalEndCapDirectionCorrectedADC;
                delete i_hHCalOtherDirectionCorrectionADC;
                delete i_file;
            }
        }
    }

    m_hHCalBarrelDirectionCorrectedADC->Sumw2();
    m_hHCalEndCapDirectionCorrectedADC->Sumw2();
    m_hHCalOtherDirectionCorrectedADC->Sumw2();

    TCanvas *pCanvas = new TCanvas("Canvas", "Canvas", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hHCalBarrelDirectionCorrectedADC->Draw();

    pCanvas->cd(2);
    m_hHCalEndCapDirectionCorrectedADC->Draw();

    pCanvas->cd(3);
    m_hHCalOtherDirectionCorrectedADC->Draw();

    TString pngOutputFilename = m_outputPath + "ADC_Distributions_and_Direction_Correction_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.png";
    TString dotCOutputFilename = m_outputPath + "ADC_Distributions_and_Direction_Correction_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.C";

    pCanvas->SaveAs(pngOutputFilename);
    pCanvas->SaveAs(dotCOutputFilename);

    m_PeakHCalBarrelDirectionCorrectedADC = m_hHCalBarrelDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalBarrelDirectionCorrectedADC));
    m_PeakHCalEndCapDirectionCorrectedADC = m_hHCalEndCapDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalEndCapDirectionCorrectedADC));
    m_PeakHCalOtherDirectionCorrectedADC = m_hHCalOtherDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalOtherDirectionCorrectedADC));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int HCalRingDigitisation::PeakFinder(const TH1F *const pTH1F)
{
    TH1F *p_TH1F = const_cast<TH1F *>(pTH1F);

    const unsigned int nBinsX(p_TH1F->GetNbinsX());

    float   previousBin;
    float   currentBin;
    float   nextBin;

    float   previousBinContent;
    float   currentBinContent;
    float   nextBinContent;

    float   highestPeak = 0.f;
    int     highestPeakBin = 0;

    for (unsigned int i = 11; i < nBinsX-10; i++)
    {
        previousBin = i-10;
        currentBin = i;
        nextBin = i+10;

        previousBinContent = p_TH1F->GetBinContent(previousBin);
        currentBinContent = p_TH1F->GetBinContent(currentBin);
        nextBinContent = p_TH1F->GetBinContent(nextBin);

        if((currentBinContent-nextBinContent) > 0 && (currentBinContent-previousBinContent) > 0 && (currentBinContent > highestPeak))
        {
            highestPeak = currentBinContent;
            highestPeakBin = i;
        }
    }
    return highestPeakBin;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], HCalRingDigitisation &hCalRingDigitisation)
{
    int c(0);

    while (((c = getopt(argc, argv, "e:i:o:f:h")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'e':
            hCalRingDigitisation.m_trueEnergy = atof(optarg);
            break;
        case 'i':
            hCalRingDigitisation.m_inputMuonRootFiles = optarg;
            break;
        case 'o':
            hCalRingDigitisation.m_outputPath = optarg;
            break;
        case 'h':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -e value  (mandatory, true energy of muons being used for calibration)" << std::endl
                      << "    -i        (mandatory, input file name(s), can include wildcards if string is in quotes)" << std::endl
                      << "    -o value  (mandatory, output path to send results to)" << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
