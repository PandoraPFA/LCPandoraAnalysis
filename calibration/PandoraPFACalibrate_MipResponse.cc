/**
 *  @file   PandoraAnalysis/calibration/PandoraPFACalibrate_MipResponse.cc
 * 
 *  @brief  Finds peak in direction corrected calo hit energy.  Used for setting MIP scale in PandoraPFA (ECalGeVToMIP, HCalGeVToMIP 
 *          and MuonGeVToMIP).
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
 *  @brief  GeVToMIP class
 */
class GeVToMIP 
{
public:
    /**
     *  @brief  Constructor
    */
    GeVToMIP();

    /**
     *  @brief  Destructor
     */
    ~GeVToMIP();

    /**
     *  @brief  Create histograms
    */
    void MakeHistograms();

// Non trivial setting on initialisation
    std::string     m_inputMuonRootFiles;                       ///< Input root files 
    float           m_trueEnergy;                               ///< True energy (opposed to kinetic) of particle being simulated
    std::string     m_outputPath;                               ///< Output path to send results
    float           m_PeakECalDirectionCorrectedCaloHitEnergy;  ///< Peak position in m_hECalDirectionCorrectedCaloHitEnergy
    float           m_PeakHCalDirectionCorrectedCaloHitEnergy;  ///< Peak position in m_hHCalDirectionCorrectedCaloHitEnergy
    float           m_PeakMuonDirectionCorrectedCaloHitEnergy;  ///< Peak position in m_hMuonDirectionCorrectedCaloHitEnergy

private:
    /**
     *  @brief  Histogram peak finder
    */
    int PeakFinder(const TH1F *const pTH1F);

    TH1F           *m_hECalDirectionCorrectedCaloHitEnergy;     ///< Histogram of direction corrections applied to calo hit energy in ECal
    TH1F           *m_hHCalDirectionCorrectedCaloHitEnergy;     ///< Histogram of direction corrections applied to calo hit energy in HCal
    TH1F           *m_hMuonDirectionCorrectedCaloHitEnergy;     ///< Histogram of direction corrections applied to calo hit energy in Muon
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

bool ParseCommandLine(int argc, char *argv[], GeVToMIP &geVToMIP);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        GeVToMIP geVToMIP;

        if (!ParseCommandLine(argc, argv, geVToMIP))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        geVToMIP.MakeHistograms();

        std::string dataFileName(geVToMIP.m_outputPath + "Calibration.txt");
        ofstream    data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "_____________________________________________________________________________________" << std::endl;
        std::cout << "_____________________________________________________________________________________" << std::endl;

        data_file << "GeV To MIP Calibration for the ECal                : " << std::endl << std::endl;
        data_file << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " << std::endl;
        data_file << "The direction corrected calo hit energy has a peak : " << std::endl;
        data_file << "in the ECal at:                                    : " << std::endl;
        data_file << "ECal GeV To MIP Peak                               : " << geVToMIP.m_PeakECalDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        data_file << "ECalGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakECalDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

        data_file << "GeV To MIP Calibration for the HCal                : " << std::endl << std::endl;
        data_file << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " << std::endl;
        data_file << "The direction corrected calo hit energy has a peak : " << std::endl;
        data_file << "in the HCal at:                                    : " << std::endl;
        data_file << "HCal GeV To MIP Peak                               : " << geVToMIP.m_PeakHCalDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        data_file << "HCalGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakHCalDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

        data_file << "GeV To MIP Calibration for the Muon Chamber        : " << std::endl << std::endl;
        data_file << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " << std::endl;
        data_file << "The direction corrected calo hit energy has a peak : " << std::endl;
        data_file << "in the Muon Chamber at:                            : " << std::endl;
        data_file << "Muon GeV To MIP Peak                               : " << geVToMIP.m_PeakMuonDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        data_file << "MuonGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakMuonDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

        std::cout << "GeV To MIP Calibration for the ECal                : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " <<std::endl;
        std::cout << "The direction corrected calo hit energy has a peak : " << std::endl;
        std::cout << "in the ECal at:                                    : " << std::endl;
        std::cout << "ECal GeV To MIP Peak                               : " << geVToMIP.m_PeakECalDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        std::cout << "ECalGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakECalDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

        std::cout << "GeV To MIP Calibration for the HCal                : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " << std::endl;
        std::cout << "The direction corrected calo hit energy has a peak : " << std::endl;
        std::cout << "in the HCal at:                                    : " << std::endl;
        std::cout << "HCal GeV To MIP Peak                               : " << geVToMIP.m_PeakHCalDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        std::cout << "HCalGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakHCalDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

        std::cout << "GeV To MIP Calibration for the Muon Chamber        : " << std::endl << std::endl;
        std::cout << "For Muons with energy                              : " << geVToMIP.m_trueEnergy << " : " << std::endl;
        std::cout << "The direction corrected calo hit energy has a peak : " << std::endl;
        std::cout << "in the Muon Chamber at:                            : " << std::endl;
        std::cout << "Muon GeV To MIP Peak                               : " << geVToMIP.m_PeakMuonDirectionCorrectedCaloHitEnergy << " : " <<std::endl;
        std::cout << "MuonGeVToMIP                                       : " << 1.f/geVToMIP.m_PeakMuonDirectionCorrectedCaloHitEnergy << " : " << std::endl << std::endl;

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

GeVToMIP::GeVToMIP() :
    m_inputMuonRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_outputPath(""),
    m_PeakECalDirectionCorrectedCaloHitEnergy(std::numeric_limits<float>::max()),
    m_PeakHCalDirectionCorrectedCaloHitEnergy(std::numeric_limits<float>::max()),
    m_PeakMuonDirectionCorrectedCaloHitEnergy(std::numeric_limits<float>::max()),
    m_hECalDirectionCorrectedCaloHitEnergy(NULL),
    m_hHCalDirectionCorrectedCaloHitEnergy(NULL),
    m_hMuonDirectionCorrectedCaloHitEnergy(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

GeVToMIP::~GeVToMIP()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GeVToMIP::MakeHistograms()
{
    m_hECalDirectionCorrectedCaloHitEnergy = new TH1F("ECalDirectionCorrectedCaloHitEnergy", "ECalDirectionCorrectedCaloHitEnergy : 1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 0.1);
    m_hECalDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in ECal");
    m_hECalDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");

    m_hHCalDirectionCorrectedCaloHitEnergy = new TH1F("HCalDirectionCorrectedCaloHitEnergy", "HCalDirectionCorrectedCaloHitEnergy : 1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 0.1);
    m_hHCalDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in HCal");
    m_hHCalDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");

    m_hMuonDirectionCorrectedCaloHitEnergy = new TH1F("MuonDirectionCorrectedCaloHitEnergy", "MuonDirectionCorrectedCaloHitEnergy : 1==nPfoTargetsTotal && 1==nPfoTargetsTracks", 500, 0., 1.0);
    m_hMuonDirectionCorrectedCaloHitEnergy->GetXaxis()->SetTitle("Direction Corrected Calo Hit Energy in Muon Chamber");
    m_hMuonDirectionCorrectedCaloHitEnergy->GetYaxis()->SetTitle("Entries");

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
                TH1F *i_hECalDirectionCorrectedCaloHitEnergy = (TH1F*) i_file->Get("ECalDirectionCorrectedCaloHitEnergy");
                TH1F *i_hHCalDirectionCorrectedCaloHitEnergy = (TH1F*) i_file->Get("HCalDirectionCorrectedCaloHitEnergy");
                TH1F *i_hMuonDirectionCorrectedCaloHitEnergy = (TH1F*) i_file->Get("MuonDirectionCorrectedCaloHitEnergy");

                if (i_hECalDirectionCorrectedCaloHitEnergy!=NULL)
                    m_hECalDirectionCorrectedCaloHitEnergy->Add(i_hECalDirectionCorrectedCaloHitEnergy,1.0);

                if (i_hHCalDirectionCorrectedCaloHitEnergy!=NULL)
                    m_hHCalDirectionCorrectedCaloHitEnergy->Add(i_hHCalDirectionCorrectedCaloHitEnergy,1.0);

                if (i_hMuonDirectionCorrectedCaloHitEnergy!=NULL)
                    m_hMuonDirectionCorrectedCaloHitEnergy->Add(i_hMuonDirectionCorrectedCaloHitEnergy,1.0);

                delete i_hECalDirectionCorrectedCaloHitEnergy;
                delete i_hHCalDirectionCorrectedCaloHitEnergy;
                delete i_hMuonDirectionCorrectedCaloHitEnergy;
                delete i_file;
            }
        }
    }

    m_hECalDirectionCorrectedCaloHitEnergy->Sumw2();
    m_hHCalDirectionCorrectedCaloHitEnergy->Sumw2();
    m_hMuonDirectionCorrectedCaloHitEnergy->Sumw2();

    TCanvas *pCanvas = new TCanvas("Canvas", "Canvas", 5000, 5000);
    pCanvas->Divide(1,3);

    pCanvas->cd(1);
    m_hECalDirectionCorrectedCaloHitEnergy->Draw();

    pCanvas->cd(2);
    m_hHCalDirectionCorrectedCaloHitEnergy->Draw();

    pCanvas->cd(3);
    m_hMuonDirectionCorrectedCaloHitEnergy->Draw();

    TCanvas *pCanvas1 = new TCanvas("Canvas1","Canvas1",5000,5000);
    pCanvas1->cd();
    m_hECalDirectionCorrectedCaloHitEnergy->Draw();
    int PeakECal = m_hECalDirectionCorrectedCaloHitEnergy->GetMaximumBin();
    TString pngOutputFilename1 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_ECal.png";
    TString dotCOutputFilename1 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_ECal.C";
    pCanvas1->SaveAs(pngOutputFilename1);
    pCanvas1->SaveAs(dotCOutputFilename1);

    TCanvas *pCanvas2 = new TCanvas("Canvas2","Canvas2",5000,5000);
    pCanvas2->cd();
    m_hHCalDirectionCorrectedCaloHitEnergy->Draw();
    int PeakHCal = m_hHCalDirectionCorrectedCaloHitEnergy->GetMaximumBin();
    TString pngOutputFilename2 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_HCal.png";
    TString dotCOutputFilename2 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_HCal.C";
    pCanvas2->SaveAs(pngOutputFilename2);
    pCanvas2->SaveAs(dotCOutputFilename2);

    TCanvas *pCanvas3 = new TCanvas("Canvas3","Canvas3",5000,5000);
    pCanvas3->cd();
    m_hMuonDirectionCorrectedCaloHitEnergy->Draw();
    int PeakMuon = m_hMuonDirectionCorrectedCaloHitEnergy->GetMaximumBin();
    TString pngOutputFilename3 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_Muon_Chamber.png";
    TString dotCOutputFilename3 = m_outputPath + "GeVToMIP_Calibration_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Muons_Muon_Chamber.C";
    pCanvas3->SaveAs(pngOutputFilename3);
    pCanvas3->SaveAs(dotCOutputFilename3);

    m_PeakECalDirectionCorrectedCaloHitEnergy = m_hECalDirectionCorrectedCaloHitEnergy->GetXaxis()->GetBinCenter(PeakECal);
    m_PeakHCalDirectionCorrectedCaloHitEnergy = m_hHCalDirectionCorrectedCaloHitEnergy->GetXaxis()->GetBinCenter(PeakHCal);
    m_PeakMuonDirectionCorrectedCaloHitEnergy = m_hMuonDirectionCorrectedCaloHitEnergy->GetXaxis()->GetBinCenter(PeakMuon);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int GeVToMIP::PeakFinder(const TH1F *const pTH1F)
{
    TH1F *p_TH1F = const_cast<TH1F *>(pTH1F);
    const unsigned int nBinsX(p_TH1F->GetNbinsX());

    float previousBin;
    float currentBin;
    float nextBin;

    float previousBinContent;
    float currentBinContent;
    float nextBinContent;

    float highestPeak = 0.f;
    int   highestPeakBin = 0;

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

bool ParseCommandLine(int argc, char *argv[], GeVToMIP &geVToMIP)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            geVToMIP.m_inputMuonRootFiles = optarg;
            break;
        case 'b':
            geVToMIP.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            geVToMIP.m_outputPath = optarg;
            break;
        case 'd':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)   " << std::endl
                      << "    -b value  (mandatory, true energy of muons being used for calibration)                    " << std::endl
                      << "    -c value  (mandatory, output path to send results to)                                     " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
