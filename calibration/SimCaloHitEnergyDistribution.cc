/**
 *  @file   PandoraAnalysis/calibration/SimCaloHitEnergyDistribution.cc
 * 
 *  @brief  Mean of distribution of direction corrections applied to sim calorimeter hits in HCal Barrel, EndCap and Other 
 *          for HCal ring digitisation.
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TSystemDirectory.h"

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  SimCaloHitEnergyDistribution class
 */
class SimCaloHitEnergyDistribution 
{
public:
    /**
     *  @brief  Constructor
    */
    SimCaloHitEnergyDistribution();

    /**
     *  @brief  Destructor
     */
    ~SimCaloHitEnergyDistribution();

    /**
     *  @brief  Create histograms
    */
    void MakeHistograms();

// Inputs Set By Parsing Command Line
    std::string     m_inputMuonRootFiles;                   ///< Input root files - Muons
    float           m_trueEnergy;                           ///< True energy (opposed to kinetic) of particle being simulated
    std::string     m_outputPath;                           ///< Output path to send results

// Outputs
    float           m_PeakHCalBarrelDirectionCorrectedSimCaloHit;  ///< Peak position in m_hHCalBarrelDirectionCorrectedSimCaloHit
    float           m_PeakHCalEndCapDirectionCorrectedSimCaloHit;  ///< Peak position in m_hHCalEndCapDirectionCorrectedSimCaloHit
    float           m_PeakHCalOtherDirectionCorrectedSimCaloHit;   ///< Peak position in m_hHCalOtherDirectionCorrectedSimCaloHit
    float           m_PeakECalDirectionCorrectedSimCaloHit;        ///< Peak position in m_hECalDirectionCorrectedSimCaloHit

private:
    /**
     *  @brief  Histogram peak finder
    */
    int PeakFinder(const TH1F *const pTH1F);

    TH1F           *m_hHCalBarrelDirectionCorrectedSimCaloHit;     ///< Histogram of direction corrected SimCalorimeterHit energy (SimCaloHit) in HCal Barrel
    TH1F           *m_hHCalEndCapDirectionCorrectedSimCaloHit;     ///< Histogram of direction corrected SimCalorimeterHit energy (SimCaloHit) in HCal EndCap
    TH1F           *m_hHCalOtherDirectionCorrectedSimCaloHit;      ///< Histogram of direction corrected SimCalorimeterHit energy (SimCaloHit) in HCal Other
    TH1F           *m_hECalDirectionCorrectedSimCaloHit;           ///< Histogram of direction corrected SimCalorimeterHit energy (SimCaloHit) in ECal
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

bool ParseCommandLine(int argc, char *argv[], SimCaloHitEnergyDistribution &simCaloHitEnergyDistribution);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        SimCaloHitEnergyDistribution simCaloHitEnergyDistribution;

        if (!ParseCommandLine(argc, argv, simCaloHitEnergyDistribution))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        simCaloHitEnergyDistribution.MakeHistograms();

        std::string dataFileName(simCaloHitEnergyDistribution.m_outputPath + "Calibration.txt");
        ofstream    data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "_____________________________________________________________________________________" << std::endl;
        std::cout << "_____________________________________________________________________________________" << std::endl;

        data_file << "Digitisation of the HCal Ring                      : " << std::endl << std::endl;
        data_file << "For Muons with energy                              : " << simCaloHitEnergyDistribution.m_trueEnergy << " : " <<std::endl;
        data_file << "The MIP peaks for each part of the HCal were found : " << std::endl;
        data_file << "to occur at the following SimCaloHit values:              : " << std::endl;
        data_file << "HCal Barrel MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalBarrelDirectionCorrectedSimCaloHit << " : " <<std::endl;
        data_file << "HCal EndCap MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalEndCapDirectionCorrectedSimCaloHit << " : " <<std::endl;
        data_file << "HCal Ring MIP Peak:                                : " << simCaloHitEnergyDistribution.m_PeakHCalOtherDirectionCorrectedSimCaloHit << " : " <<std::endl;
        data_file << "ECal MIP Peak:                                     : " << simCaloHitEnergyDistribution.m_PeakECalDirectionCorrectedSimCaloHit << " : " <<std::endl;

        std::cout << "Digitisation of the HCal Ring                      : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << simCaloHitEnergyDistribution.m_trueEnergy << " : " <<std::endl;
        std::cout << "The MIP peaks for each part of the HCal were found : " << std::endl;
        std::cout << "to occur at the following SimCaloHit values:              : " << std::endl;
        std::cout << "HCal Barrel MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalBarrelDirectionCorrectedSimCaloHit << " : " <<std::endl;
        std::cout << "HCal EndCap MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalEndCapDirectionCorrectedSimCaloHit << " : " <<std::endl;
        std::cout << "HCal Ring MIP Peak:                                : " << simCaloHitEnergyDistribution.m_PeakHCalOtherDirectionCorrectedSimCaloHit << " : " <<std::endl;
        std::cout << "ECal MIP Peak:                                     : " << simCaloHitEnergyDistribution.m_PeakECalDirectionCorrectedSimCaloHit << " : " <<std::endl;

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

SimCaloHitEnergyDistribution::SimCaloHitEnergyDistribution() :
    m_inputMuonRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_outputPath(""),
    m_PeakHCalBarrelDirectionCorrectedSimCaloHit(std::numeric_limits<float>::max()),
    m_PeakHCalEndCapDirectionCorrectedSimCaloHit(std::numeric_limits<float>::max()),
    m_PeakHCalOtherDirectionCorrectedSimCaloHit(std::numeric_limits<float>::max()),
    m_PeakECalDirectionCorrectedSimCaloHit(std::numeric_limits<float>::max()),
    m_hHCalBarrelDirectionCorrectedSimCaloHit(NULL),
    m_hHCalEndCapDirectionCorrectedSimCaloHit(NULL),
    m_hHCalOtherDirectionCorrectedSimCaloHit(NULL),
    m_hECalDirectionCorrectedSimCaloHit(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimCaloHitEnergyDistribution::~SimCaloHitEnergyDistribution()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimCaloHitEnergyDistribution::MakeHistograms()
{
    m_hHCalBarrelDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitBarrel", "Distribution of Direction Corrected SimCaloHits in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalBarrelDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalBarrelDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitEndCap", "Distribution of Direction Corrected SimCaloHits in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hHCalEndCapDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalEndCapDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectedSimCaloHit = new TH1F("HCalDirectionCorrectedSimCaloHitOther", "Distribution of Direction Corrected SimCaloHits in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.005);
    m_hHCalOtherDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hHCalOtherDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalDirectionCorrectedSimCaloHit = new TH1F("ECalDirectionCorrectedSimCaloHit", "Distribution of Direction Corrected SimCaloHits in the ECal (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hECalDirectionCorrectedSimCaloHit->GetXaxis()->SetTitle("Direction Corrected SimCaloHit Measurement");
    m_hECalDirectionCorrectedSimCaloHit->GetYaxis()->SetTitle("Entries");

    unsigned int found_slash = m_inputMuonRootFiles.find_last_of("/");
    std::string path = m_inputMuonRootFiles.substr(0,found_slash);
    std::string file = m_inputMuonRootFiles.substr(found_slash+1);

    unsigned int found_star = file.find_last_of("*");
    std::string file_prefix = file.substr(0,found_star);
    std::string file_suffix = file.substr(found_star+1);

    TSystemDirectory dir(path.c_str(), path.c_str());
    TList *pTList = dir.GetListOfFiles();

    if (pTList)
    {
        TSystemFile *pTSystemFile;
        TString fname;
        TIter next(pTList);
        while ((pTSystemFile = (TSystemFile*)next()))
        {
            fname = pTSystemFile->GetName();

            if (!pTSystemFile->IsDirectory() && fname.EndsWith(file_suffix.c_str()) && fname.BeginsWith(file_prefix.c_str()))
            {
                TString filename(path + "/" + fname);
                TFile *pTFile = new TFile(filename);
                TH1F *pTH1FHCalBarrelDirectionCorrectedSimCaloHit = (TH1F*) pTFile->Get("HCalDirectionCorrectedSimCaloHitBarrel");
                TH1F *pTH1FHCalEndCapDirectionCorrectedSimCaloHit = (TH1F*) pTFile->Get("HCalDirectionCorrectedSimCaloHitEndCap");
                TH1F *pTH1FHCalOtherDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("HCalDirectionCorrectedSimCaloHitOther");
                TH1F *pTH1FECalDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("ECalDirectionCorrectedSimCaloHit");

                if (pTH1FHCalBarrelDirectionCorrectedSimCaloHit!=NULL)
                    m_hHCalBarrelDirectionCorrectedSimCaloHit->Add(pTH1FHCalBarrelDirectionCorrectedSimCaloHit,1.0);

                if (pTH1FHCalEndCapDirectionCorrectedSimCaloHit!=NULL)
                    m_hHCalEndCapDirectionCorrectedSimCaloHit->Add(pTH1FHCalEndCapDirectionCorrectedSimCaloHit,1.0);

                if (pTH1FHCalOtherDirectionCorrectionSimCaloHit!=NULL)
                    m_hHCalOtherDirectionCorrectedSimCaloHit->Add(pTH1FHCalOtherDirectionCorrectionSimCaloHit,1.0);

                if (pTH1FECalDirectionCorrectionSimCaloHit!=NULL)
                    m_hECalDirectionCorrectedSimCaloHit->Add(pTH1FECalDirectionCorrectionSimCaloHit,1.0);

                delete pTH1FHCalBarrelDirectionCorrectedSimCaloHit;
                delete pTH1FHCalEndCapDirectionCorrectedSimCaloHit;
                delete pTH1FHCalOtherDirectionCorrectionSimCaloHit;
                delete pTH1FECalDirectionCorrectionSimCaloHit;
                delete pTFile;
            }
        }
    }

    m_hHCalBarrelDirectionCorrectedSimCaloHit->Sumw2();
    m_hHCalEndCapDirectionCorrectedSimCaloHit->Sumw2();
    m_hHCalOtherDirectionCorrectedSimCaloHit->Sumw2();
    m_hECalDirectionCorrectedSimCaloHit->Sumw2();

    TCanvas *pCanvas = new TCanvas("Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal", "Direction Corrected SimCalorimeterHit Energy Distribution HCal", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hHCalBarrelDirectionCorrectedSimCaloHit->Draw();

    pCanvas->cd(2);
    m_hHCalEndCapDirectionCorrectedSimCaloHit->Draw();

    pCanvas->cd(3);
    m_hHCalOtherDirectionCorrectedSimCaloHit->Draw();

    TString pngOutputFilename = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.png";
    TString dotCOutputFilename = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.C";

    pCanvas->SaveAs(pngOutputFilename);
    pCanvas->SaveAs(dotCOutputFilename);

    TCanvas *pCanvas2 = new TCanvas("Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal", "Direction Corrected SimCalorimeterHit Energy Distribution ECal", 5000, 5000);
    m_hECalDirectionCorrectedSimCaloHit->Draw();

    TString pngOutputFilename2 = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.png";
    TString dotCOutputFilename2 = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.C";

    pCanvas2->SaveAs(pngOutputFilename2);
    pCanvas2->SaveAs(dotCOutputFilename2);

    m_PeakHCalBarrelDirectionCorrectedSimCaloHit = m_hHCalBarrelDirectionCorrectedSimCaloHit->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalBarrelDirectionCorrectedSimCaloHit));
    m_PeakHCalEndCapDirectionCorrectedSimCaloHit = m_hHCalEndCapDirectionCorrectedSimCaloHit->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalEndCapDirectionCorrectedSimCaloHit));
    m_PeakHCalOtherDirectionCorrectedSimCaloHit = m_hHCalOtherDirectionCorrectedSimCaloHit->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalOtherDirectionCorrectedSimCaloHit));
    m_PeakECalDirectionCorrectedSimCaloHit = m_hECalDirectionCorrectedSimCaloHit->GetXaxis()->GetBinCenter(PeakFinder(m_hECalDirectionCorrectedSimCaloHit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int SimCaloHitEnergyDistribution::PeakFinder(const TH1F *const pTH1F)
{
    TH1F *p_TH1F = const_cast<TH1F *>(pTH1F);

    const unsigned int nBinsX(p_TH1F->GetNbinsX());

    float previousBin, currentBin, nextBin, previousBinContent, currentBinContent, nextBinContent, highestPeak = 0.f;
    int highestPeakBin = 0;

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

bool ParseCommandLine(int argc, char *argv[], SimCaloHitEnergyDistribution &simCaloHitEnergyDistribution)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            simCaloHitEnergyDistribution.m_inputMuonRootFiles = optarg;
            break;
        case 'b':
            simCaloHitEnergyDistribution.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            simCaloHitEnergyDistribution.m_outputPath = optarg;
            break;
        case 'd':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input muon file name(s), can include wildcards if string is in quotes)  " << std::endl
                      << "    -b value  (mandatory, true energy of muons being used for calibration)                        " << std::endl
                      << "    -c value  (mandatory, output path to send results to)                                         " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
