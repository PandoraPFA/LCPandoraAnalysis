/**
 *  @file   PandoraAnalysis/calibration/ECalDigitisation_DirectionCorrectionDistribution.cc
 *
 *  @brief  Generate average direction correction applied to photon events in ECal Barrel, EndCap and Other.  Used for setting
 *          digitisation constant CalibrECalOther.
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
 *  @brief  ECalDigitisationDCDistribution class
 */
class ECalDigitisationDCDistribution
{
public:
    /**
     *  @brief  Constructor
    */
    ECalDigitisationDCDistribution();

    /**
     *  @brief  Destructor
     */
    ~ECalDigitisationDCDistribution();

    /**
     *  @brief  Create histograms
    */
    void MakeHistograms();

// Inputs Set By Parsing Command Line
    std::string     m_inputPhotonRootFiles;                  ///< Input root files - Photon
    float           m_trueEnergy;                           ///< True energy (opposed to kinetic) of particle being simulated
    std::string     m_outputPath;                           ///< Output path to send results

// Outputs
    TH1F           *m_hECalBarrelDirectionCorrectionSimCaloHit;    ///< Histogram of direction corrections applied to SimCalorimeterHits in ECal Barrel
    TH1F           *m_hECalEndCapDirectionCorrectionSimCaloHit;    ///< Histogram of direction corrections applied to SimCalorimeterHits in ECal EndCap
    TH1F           *m_hECalOtherDirectionCorrectionSimCaloHit;     ///< Histogram of direction corrections applied to SimCalorimeterHits in ECal Other
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

bool ParseCommandLine(int argc, char *argv[], ECalDigitisationDCDistribution &eCalDigitisationDCDistribution);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        ECalDigitisationDCDistribution eCalDigitisationDCDistribution;

        if (!ParseCommandLine(argc, argv, eCalDigitisationDCDistribution))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        eCalDigitisationDCDistribution.MakeHistograms();

        std::string dataFileName(eCalDigitisationDCDistribution.m_outputPath + "Calibration.txt");
        std::ofstream data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "The average direction correction applied to Photon  : " << std::endl;
        data_file << "events with energy:                                : " << eCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        data_file << "in the ECal EndCap is given by:                    : " << std::endl;
        data_file << "Mean Direction Correction ECalEndCap:              : " << eCalDigitisationDCDistribution.m_hECalEndCapDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to Photon  : " << std::endl;
        std::cout << "events with energy:                                : " << eCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the ECal EndCap is given by:                    : " << std::endl;
        std::cout << "Mean Direction Correction ECalEndCap:              : " << eCalDigitisationDCDistribution.m_hECalEndCapDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        data_file << "The average direction correction applied to Photon  : " << std::endl;
        data_file << "events with energy:                                : " << eCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        data_file << "in the ECalOther (i.e. ECal Ring) is given by:     : " << std::endl;
        data_file << "Mean Direction Correction ECalOther:               : " << eCalDigitisationDCDistribution.m_hECalOtherDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

        std::cout << "The average direction correction applied to Photon  : " << std::endl;
        std::cout << "events with energy:                                : " << eCalDigitisationDCDistribution.m_trueEnergy << " : /GeV" <<std::endl;
        std::cout << "in the ECalOther (i.e. ECal Ring) is given by:     : " << std::endl;
        std::cout << "Mean Direction Correction ECalOther:               : " << eCalDigitisationDCDistribution.m_hECalOtherDirectionCorrectionSimCaloHit->GetMean() << " : " <<std::endl<<std::endl;

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

ECalDigitisationDCDistribution::ECalDigitisationDCDistribution() :
    m_inputPhotonRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_outputPath(""),
    m_hECalBarrelDirectionCorrectionSimCaloHit(NULL),
    m_hECalEndCapDirectionCorrectionSimCaloHit(NULL),
    m_hECalOtherDirectionCorrectionSimCaloHit(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ECalDigitisationDCDistribution::~ECalDigitisationDCDistribution()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ECalDigitisationDCDistribution::MakeHistograms()
{
    m_hECalBarrelDirectionCorrectionSimCaloHit = new TH1F("ECalBarrelDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal Barrel (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in ECal)", 200, 0., 1.0);
    m_hECalBarrelDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hECalBarrelDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalEndCapDirectionCorrectionSimCaloHit = new TH1F("ECalEndCapDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal EndCap (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in ECal)", 200, 0., 1.0);
    m_hECalEndCapDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hECalEndCapDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    m_hECalOtherDirectionCorrectionSimCaloHit = new TH1F("ECalOtherDirectionCorrectionSimCaloHit", "Distribution of Direction Corrections for SimCaloHits in the ECal Other (1==nPfoTargetsTotal && 1==nPfoTargetsNeutralHadrons && Contained in ECal)", 200, 0., 1.0);
    m_hECalOtherDirectionCorrectionSimCaloHit->GetXaxis()->SetTitle("Sim Calo Hit Direction Correction");
    m_hECalOtherDirectionCorrectionSimCaloHit->GetYaxis()->SetTitle("Entries");

    unsigned int found_slash = m_inputPhotonRootFiles.find_last_of("/");
    std::string path = m_inputPhotonRootFiles.substr(0,found_slash);
    std::string file = m_inputPhotonRootFiles.substr(found_slash+1);

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
                TH1F *pTH1FECalBarrelDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("ECalBarrelDirectionCorrectionSimCaloHit");
                TH1F *pTH1FECalEndCapDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("ECalEndCapDirectionCorrectionSimCaloHit");
                TH1F *pTH1FECalOtherDirectionCorrectionSimCaloHit = (TH1F*) pTFile->Get("ECalOtherDirectionCorrectionSimCaloHit");

                if (pTH1FECalBarrelDirectionCorrectionSimCaloHit!=NULL)
                    m_hECalBarrelDirectionCorrectionSimCaloHit->Add(pTH1FECalBarrelDirectionCorrectionSimCaloHit,1.0);

                if (pTH1FECalEndCapDirectionCorrectionSimCaloHit!=NULL)
                    m_hECalEndCapDirectionCorrectionSimCaloHit->Add(pTH1FECalEndCapDirectionCorrectionSimCaloHit,1.0);

                if (pTH1FECalOtherDirectionCorrectionSimCaloHit!=NULL)
                    m_hECalOtherDirectionCorrectionSimCaloHit->Add(pTH1FECalOtherDirectionCorrectionSimCaloHit,1.0);

                delete pTH1FECalBarrelDirectionCorrectionSimCaloHit;
                delete pTH1FECalEndCapDirectionCorrectionSimCaloHit;
                delete pTH1FECalOtherDirectionCorrectionSimCaloHit;
                delete pTFile;
            }
        }
    }

    m_hECalBarrelDirectionCorrectionSimCaloHit->Sumw2();
    m_hECalEndCapDirectionCorrectionSimCaloHit->Sumw2();
    m_hECalOtherDirectionCorrectionSimCaloHit->Sumw2();

    TCanvas *pCanvas = new TCanvas("Canvas", "Canvas", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hECalBarrelDirectionCorrectionSimCaloHit->Draw();

    pCanvas->cd(2);
    m_hECalEndCapDirectionCorrectionSimCaloHit->Draw();

    pCanvas->cd(3);
    m_hECalOtherDirectionCorrectionSimCaloHit->Draw();

    TString pngOutputFilename = m_outputPath + "Direction_Correction_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Photon.png";
    TString dotCOutputFilename = m_outputPath + "Direction_Correction_Distribution_ECal_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Photon.C";

    pCanvas->SaveAs(pngOutputFilename);
    pCanvas->SaveAs(dotCOutputFilename);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ECalDigitisationDCDistribution &eCalDigitisationDCDistribution)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            eCalDigitisationDCDistribution.m_inputPhotonRootFiles = optarg;
            break;
        case 'b':
            eCalDigitisationDCDistribution.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            eCalDigitisationDCDistribution.m_outputPath = optarg;
            break;
        case 'd':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)   " << std::endl
                      << "    -b value  (mandatory, true energy of photon being used for calibration)                    " << std::endl
                      << "    -c        (mandatory, output path to send results to)                                     " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
