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
    float           m_PeakHCalBarrelDirectionCorrectedADC;  ///< Peak position in m_hHCalBarrelDirectionCorrectedADC
    float           m_PeakHCalEndCapDirectionCorrectedADC;  ///< Peak position in m_hHCalEndCapDirectionCorrectedADC
    float           m_PeakHCalOtherDirectionCorrectedADC;   ///< Peak position in m_hHCalOtherDirectionCorrectedADC
    float           m_PeakECalDirectionCorrectedADC;        ///< Peak position in m_hECalDirectionCorrectedADC

private:
    /**
     *  @brief  Histogram peak finder
    */
    int PeakFinder(const TH1F *const pTH1F);

    TH1F           *m_hHCalBarrelDirectionCorrectedADC;     ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal Barrel
    TH1F           *m_hHCalEndCapDirectionCorrectedADC;     ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal EndCap
    TH1F           *m_hHCalOtherDirectionCorrectedADC;      ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in HCal Other
    TH1F           *m_hECalDirectionCorrectedADC;           ///< Histogram of direction corrected SimCalorimeterHit energy (ADC) in ECal
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
        data_file << "to occur at the following ADC values:              : " << std::endl;
        data_file << "HCal Barrel MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalBarrelDirectionCorrectedADC << " : " <<std::endl;
        data_file << "HCal EndCap MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalEndCapDirectionCorrectedADC << " : " <<std::endl;
        data_file << "HCal Ring MIP Peak:                                : " << simCaloHitEnergyDistribution.m_PeakHCalOtherDirectionCorrectedADC << " : " <<std::endl;
        data_file << "ECal MIP Peak:                                     : " << simCaloHitEnergyDistribution.m_PeakECalDirectionCorrectedADC << " : " <<std::endl;

        std::cout << "Digitisation of the HCal Ring                      : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << simCaloHitEnergyDistribution.m_trueEnergy << " : " <<std::endl;
        std::cout << "The MIP peaks for each part of the HCal were found : " << std::endl;
        std::cout << "to occur at the following ADC values:              : " << std::endl;
        std::cout << "HCal Barrel MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalBarrelDirectionCorrectedADC << " : " <<std::endl;
        std::cout << "HCal EndCap MIP Peak:                              : " << simCaloHitEnergyDistribution.m_PeakHCalEndCapDirectionCorrectedADC << " : " <<std::endl;
        std::cout << "HCal Ring MIP Peak:                                : " << simCaloHitEnergyDistribution.m_PeakHCalOtherDirectionCorrectedADC << " : " <<std::endl;
        std::cout << "ECal MIP Peak:                                     : " << simCaloHitEnergyDistribution.m_PeakECalDirectionCorrectedADC << " : " <<std::endl;

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
    m_PeakHCalBarrelDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_PeakHCalEndCapDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_PeakHCalOtherDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_PeakECalDirectionCorrectedADC(std::numeric_limits<float>::max()),
    m_hHCalBarrelDirectionCorrectedADC(NULL),
    m_hHCalEndCapDirectionCorrectedADC(NULL),
    m_hHCalOtherDirectionCorrectedADC(NULL),
    m_hECalDirectionCorrectedADC(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SimCaloHitEnergyDistribution::~SimCaloHitEnergyDistribution()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimCaloHitEnergyDistribution::MakeHistograms()
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

    m_hECalDirectionCorrectedADC = new TH1F("ECalDirectionCorrectedADC", "Distribution of Direction Corrected ADCs in the ECal (1==nPfoTargetsTotal && 1==nPfoTargetsTracks)", 200, 0., 0.001);
    m_hECalDirectionCorrectedADC->GetXaxis()->SetTitle("Direction Corrected ADC Measurement");
    m_hECalDirectionCorrectedADC->GetYaxis()->SetTitle("Entries");

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
                TH1F *pTH1FHCalBarrelDirectionCorrectedADC = (TH1F*) pTFile->Get("HCalDirectionCorrectedADCBarrel");
                TH1F *pTH1FHCalEndCapDirectionCorrectedADC = (TH1F*) pTFile->Get("HCalDirectionCorrectedADCEndCap");
                TH1F *pTH1FHCalOtherDirectionCorrectionADC = (TH1F*) pTFile->Get("HCalDirectionCorrectedADCOther");
                TH1F *pTH1FECalDirectionCorrectionADC = (TH1F*) pTFile->Get("ECalDirectionCorrectedADC");

                if (pTH1FHCalBarrelDirectionCorrectedADC!=NULL)
                    m_hHCalBarrelDirectionCorrectedADC->Add(pTH1FHCalBarrelDirectionCorrectedADC,1.0);

                if (pTH1FHCalEndCapDirectionCorrectedADC!=NULL)
                    m_hHCalEndCapDirectionCorrectedADC->Add(pTH1FHCalEndCapDirectionCorrectedADC,1.0);

                if (pTH1FHCalOtherDirectionCorrectionADC!=NULL)
                    m_hHCalOtherDirectionCorrectedADC->Add(pTH1FHCalOtherDirectionCorrectionADC,1.0);

                if (pTH1FECalDirectionCorrectionADC!=NULL)
                    m_hECalDirectionCorrectedADC->Add(pTH1FECalDirectionCorrectionADC,1.0);

                delete pTH1FHCalBarrelDirectionCorrectedADC;
                delete pTH1FHCalEndCapDirectionCorrectedADC;
                delete pTH1FHCalOtherDirectionCorrectionADC;
                delete pTH1FECalDirectionCorrectionADC;
                delete pTFile;
            }
        }
    }

    m_hHCalBarrelDirectionCorrectedADC->Sumw2();
    m_hHCalEndCapDirectionCorrectedADC->Sumw2();
    m_hHCalOtherDirectionCorrectedADC->Sumw2();
    m_hECalDirectionCorrectedADC->Sumw2();

    TCanvas *pCanvas = new TCanvas("Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal", "Direction Corrected SimCalorimeterHit Energy Distribution HCal", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hHCalBarrelDirectionCorrectedADC->Draw();

    pCanvas->cd(2);
    m_hHCalEndCapDirectionCorrectedADC->Draw();

    pCanvas->cd(3);
    m_hHCalOtherDirectionCorrectedADC->Draw();

    TString pngOutputFilename = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.png";
    TString dotCOutputFilename = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.C";

    pCanvas->SaveAs(pngOutputFilename);
    pCanvas->SaveAs(dotCOutputFilename);

    TCanvas *pCanvas2 = new TCanvas("Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal", "Direction Corrected SimCalorimeterHit Energy Distribution ECal", 5000, 5000);
    m_hECalDirectionCorrectedADC->Draw();

    TString pngOutputFilename2 = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.png";
    TString dotCOutputFilename2 = m_outputPath + "Direction_Corrected_SimCalorimeterHit_Energy_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons.C";

    pCanvas2->SaveAs(pngOutputFilename2);
    pCanvas2->SaveAs(dotCOutputFilename2);

    m_PeakHCalBarrelDirectionCorrectedADC = m_hHCalBarrelDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalBarrelDirectionCorrectedADC));
    m_PeakHCalEndCapDirectionCorrectedADC = m_hHCalEndCapDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalEndCapDirectionCorrectedADC));
    m_PeakHCalOtherDirectionCorrectedADC = m_hHCalOtherDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hHCalOtherDirectionCorrectedADC));
    m_PeakECalDirectionCorrectedADC = m_hECalDirectionCorrectedADC->GetXaxis()->GetBinCenter(PeakFinder(m_hECalDirectionCorrectedADC));
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
