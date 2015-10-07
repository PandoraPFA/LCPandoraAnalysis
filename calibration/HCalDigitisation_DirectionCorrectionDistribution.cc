/**
 *  @file   PandoraAnalysis/calibration/HCalDigitisation_DirectionCorrectionDistribution.cc
 * 
 *  @brief  Generate average direction correction applied to kaonL events in HCal Barrel, EndCap and Other.  Used for setting 
 *          digitisation constant CalibrHCALOther.
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
 *  @brief  HCalDigitisationDCDistribution class
 */
class HCalDigitisationDCDistribution 
{
public:
    /**
     *  @brief  Constructor
    */
    HCalDigitisationDCDistribution();

    /**
     *  @brief  Destructor
     */
    ~HCalDigitisationDCDistribution();

    /**
     *  @brief  Create histograms
    */
    void MakeHistograms();

// Inputs Set By Parsing Command Line
    std::string     m_inputKaonLRootFiles;                  ///< Input root files - KaonL
    float           m_trueEnergy;                           ///< True energy (opposed to kinetic) of particle being simulated
    std::string     m_outputPath;                           ///< Output path to send results

// Outputs
    TH1F           *m_hHCalBarrelDirectionCorrectionSimCaloHit;    ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal Barrel
    TH1F           *m_hHCalEndCapDirectionCorrectionSimCaloHit;    ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal EndCap
    TH1F           *m_hHCalOtherDirectionCorrectionSimCaloHit;     ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal Other
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

bool ParseCommandLine(int argc, char *argv[], HCalDigitisationDCDistribution &hCalDigitisationDCDistribution);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        HCalDigitisationDCDistribution hCalDigitisationDCDistribution;

        if (!ParseCommandLine(argc, argv, hCalDigitisationDCDistribution))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        hCalDigitisationDCDistribution.MakeHistograms();

        std::string dataFileName(hCalDigitisationDCDistribution.m_outputPath + "Calibration.txt");
        ofstream    data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "The average direction correction applied to KaonL  : " << std::endl;
        data_file << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        data_file << "in the HCal EndCap is given by:                    : " << std::endl;
        data_file << "Mean Direction Correction HCalEndCap:              : " << hCalDigitisationDCDistribution.m_hHCalEndCapDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to KaonL  : " << std::endl;
        std::cout << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the HCal EndCap is given by:                    : " << std::endl;
        std::cout << "Mean Direction Correction HCalEndCap:              : " << hCalDigitisationDCDistribution.m_hHCalEndCapDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        data_file << "The average direction correction applied to KaonL  : " << std::endl;
        data_file << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        data_file << "in the HCalOther (i.e. HCal Ring) is given by:     : " << std::endl;
        data_file << "Mean Direction Correction HCalOther:               : " << hCalDigitisationDCDistribution.m_hHCalOtherDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to KaonL  : " << std::endl;
        std::cout << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the HCalOther (i.e. HCal Ring) is given by:     : " << std::endl;
        std::cout << "Mean Direction Correction HCalOther:               : " << hCalDigitisationDCDistribution.m_hHCalOtherDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

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

HCalDigitisationDCDistribution::HCalDigitisationDCDistribution() :
    m_inputKaonLRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_outputPath(""),
    m_hHCalBarrelDirectionCorrectionSimCaloHit(NULL),
    m_hHCalEndCapDirectionCorrectionSimCaloHit(NULL),
    m_hHCalOtherDirectionCorrectionSimCaloHit(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HCalDigitisationDCDistribution::~HCalDigitisationDCDistribution()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisationDCDistribution::MakeHistograms()
{
    m_hHCalBarrelDirectionCorrectionSimCaloHit = new TH1F("HCalBarrelDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalBarrelDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalBarrelDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectionSimCaloHit = new TH1F("HCalEndCapDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalEndCapDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalEndCapDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectionSimCaloHit = new TH1F("HCalOtherDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalOtherDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalOtherDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    unsigned int found_slash = m_inputKaonLRootFiles.find_last_of("/");
    std::string path = m_inputKaonLRootFiles.substr(0,found_slash);
    std::string file = m_inputKaonLRootFiles.substr(found_slash+1);

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
                TH1F *pTH1FHCalBarrelDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("HCalBarrelDirectionCorrectionSimCaloHit");
                TH1F *pTH1FHCalEndCapDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("HCalEndCapDirectionCorrectionSimCaloHit");
                TH1F *pTH1FHCalOtherDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("HCalOtherDirectionCorrectionSimCaloHit");

                if (pTH1FHCalBarrelDirectionCorrectionSimCaloHit!=NULL)
                    m_hHCalBarrelDirectionCorrectionSimCaloHit->Add(pTH1FHCalBarrelDirectionCorrectionSimCaloHit,1.0);

                if (pTH1FHCalEndCapDirectionCorrectionSimCaloHit!=NULL)
                    m_hHCalEndCapDirectionCorrectionSimCaloHit->Add(pTH1FHCalEndCapDirectionCorrectionSimCaloHit,1.0);

                if (pTH1FHCalOtherDirectionCorrectionSimCaloHit!=NULL)
                    m_hHCalOtherDirectionCorrectionSimCaloHit->Add(pTH1FHCalOtherDirectionCorrectionSimCaloHit,1.0);

                delete pTH1FHCalBarrelDirectionCorrectionSimCaloHit;
                delete pTH1FHCalEndCapDirectionCorrectionSimCaloHit;
                delete pTH1FHCalOtherDirectionCorrectionSimCaloHit;
                delete pTFile;
            }
        }
    }

    m_hHCalBarrelDirectionCorrectionSimCaloHit->Sumw2();
    m_hHCalEndCapDirectionCorrectionSimCaloHit->Sumw2();
    m_hHCalOtherDirectionCorrectionSimCaloHit->Sumw2();

    TCanvas *pCanvas = new TCanvas("Canvas", "Canvas", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hHCalBarrelDirectionCorrectionSimCaloHit->Draw();

    pCanvas->cd(2);
    m_hHCalEndCapDirectionCorrectionSimCaloHit->Draw();

    pCanvas->cd(3);
    m_hHCalOtherDirectionCorrectionSimCaloHit->Draw();

    TString pngOutputFilename = m_outputPath + "Direction_Correction_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.png";
    TString dotCOutputFilename = m_outputPath + "Direction_Correction_Distribution_HCal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.C";

    pCanvas->SaveAs(pngOutputFilename);
    pCanvas->SaveAs(dotCOutputFilename);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], HCalDigitisationDCDistribution &hCalDigitisationDCDistribution)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            hCalDigitisationDCDistribution.m_inputKaonLRootFiles = optarg;
            break;
        case 'b':
            hCalDigitisationDCDistribution.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            hCalDigitisationDCDistribution.m_outputPath = optarg;
            break;
        case 'd':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)   " << std::endl
                      << "    -b value  (mandatory, true energy of kaonL being used for calibration)                    " << std::endl
                      << "    -c        (mandatory, output path to send results to)                                     " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
