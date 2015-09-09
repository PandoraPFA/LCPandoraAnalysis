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
    TH1F           *m_hHCalBarrelDirectionCorrectionADC;    ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal Barrel
    TH1F           *m_hHCalEndCapDirectionCorrectionADC;    ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal EndCap
    TH1F           *m_hHCalOtherDirectionCorrectionADC;     ///< Histogram of direction corrections applied to SimCalorimeterHits in HCal Other
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
        data_file << "Mean Direction Correction HCalEndCap:              : " << hCalDigitisationDCDistribution.m_hHCalEndCapDirectionCorrectionADC->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to KaonL  : " << std::endl;
        std::cout << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the HCal EndCap is given by:                    : " << std::endl;
        std::cout << "Mean Direction Correction HCalEndCap:              : " << hCalDigitisationDCDistribution.m_hHCalEndCapDirectionCorrectionADC->GetMean() << " : " <<std::endl<<std::endl;

        data_file << "The average direction correction applied to KaonL  : " << std::endl;
        data_file << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        data_file << "in the HCalOther (i.e. HCal Ring) is given by:     : " << std::endl;
        data_file << "Mean Direction Correction HCalOther:               : " << hCalDigitisationDCDistribution.m_hHCalOtherDirectionCorrectionADC->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to KaonL  : " << std::endl;
        std::cout << "events with energy:                                : " << hCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the HCalOther (i.e. HCal Ring) is given by:     : " << std::endl;
        std::cout << "Mean Direction Correction HCalOther:               : " << hCalDigitisationDCDistribution.m_hHCalOtherDirectionCorrectionADC->GetMean() << " : " <<std::endl<<std::endl;

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
    m_hHCalBarrelDirectionCorrectionADC(NULL),
    m_hHCalEndCapDirectionCorrectionADC(NULL),
    m_hHCalOtherDirectionCorrectionADC(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HCalDigitisationDCDistribution::~HCalDigitisationDCDistribution()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HCalDigitisationDCDistribution::MakeHistograms()
{
    m_hHCalBarrelDirectionCorrectionADC = new TH1F("HCalBarrelDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalBarrelDirectionCorrectionADC->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalBarrelDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

    m_hHCalEndCapDirectionCorrectionADC = new TH1F("HCalEndCapDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalEndCapDirectionCorrectionADC->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalEndCapDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

    m_hHCalOtherDirectionCorrectionADC = new TH1F("HCalOtherDirectionCorrectionADC", "Distribution of Direction Corrections for ADCs in the HCal Other (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in HCal)", 200, 0., 1.0);
    m_hHCalOtherDirectionCorrectionADC->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hHCalOtherDirectionCorrectionADC->GetYaxis()->SetTitle("Entries");

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
                TH1F *pTH1FHCalBarrelDirectionCorrectionADC = (TH1F*) pTFile->Get("HCalBarrelDirectionCorrectionADC");
                TH1F *pTH1FHCalEndCapDirectionCorrectionADC = (TH1F*) pTFile->Get("HCalEndCapDirectionCorrectionADC");
                TH1F *pTH1FHCalOtherDirectionCorrectionADC = (TH1F*) pTFile->Get("HCalOtherDirectionCorrectionADC");

                if (pTH1FHCalBarrelDirectionCorrectionADC!=NULL)
                    m_hHCalBarrelDirectionCorrectionADC->Add(pTH1FHCalBarrelDirectionCorrectionADC,1.0);

                if (pTH1FHCalEndCapDirectionCorrectionADC!=NULL)
                    m_hHCalEndCapDirectionCorrectionADC->Add(pTH1FHCalEndCapDirectionCorrectionADC,1.0);

                if (pTH1FHCalOtherDirectionCorrectionADC!=NULL)
                    m_hHCalOtherDirectionCorrectionADC->Add(pTH1FHCalOtherDirectionCorrectionADC,1.0);

                delete pTH1FHCalBarrelDirectionCorrectionADC;
                delete pTH1FHCalEndCapDirectionCorrectionADC;
                delete pTH1FHCalOtherDirectionCorrectionADC;
                delete pTFile;
            }
        }
    }

    m_hHCalBarrelDirectionCorrectionADC->Sumw2();
    m_hHCalEndCapDirectionCorrectionADC->Sumw2();
    m_hHCalOtherDirectionCorrectionADC->Sumw2();

    TCanvas *pCanvas = new TCanvas("Canvas", "Canvas", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hHCalBarrelDirectionCorrectionADC->Draw();

    pCanvas->cd(2);
    m_hHCalEndCapDirectionCorrectionADC->Draw();

    pCanvas->cd(3);
    m_hHCalOtherDirectionCorrectionADC->Draw();

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
