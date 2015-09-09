/**
 *  @file   PandoraAnalysis/calibration/PandoraPFACalibrate_HadronicScale_ChiSquareMethod.cc
 * 
 *  @brief  Chi square method.  Used for setting hadronic scale in PandoraPFA (ECalToHadGeVCalibration/HCalToHadGeVCalibration).
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  ChiSquaredMethod class
 */
class ChiSquaredMethod 
{
public:
    /**
     *  @brief  Constructor
    */
    ChiSquaredMethod();

    /**
     *  @brief  Destructor
     */
    ~ChiSquaredMethod();

    /**
     *  @brief  Perform the chi squared method for hadronic energy scale calibration in PandoraPFA
    */
    void Process();

// Inputs Set By Parsing Command Line
    std::string         m_inputKaonLRootFiles;          ///< Input root files - KaonL
    float               m_trueEnergy;                   ///< Total energy of calibration particle
    float               m_calibrationAccuracy;          ///< Fractional accuracy to calibrate H/ECalToHad to
    std::string         m_outputPath;                   ///< Output path to send plots to
    int                 m_numberHCalLayers;             ///< Number of layers in the HCal

// Outputs
    float               m_eCalToHadInterceptMinChi2;    ///< Best fit ECalToHad intercept
    float               m_hCalToHadInterceptMinChi2;    ///< Best fit HCalToHad intercept
    float               m_minChi2;                      ///< Best fit measure, chi squared of fit
    int                 m_nEntriesHist;                 ///< Number of entries in histogram

private:
    /**
     *  @brief  Create the histogram needed for the chi squared method
    */
    void CreateHistogram();

    /**
     *  @brief  Fill the histogram needed for the chi squared method
    */
    void FillHistogram();

    /**
     *  @brief  Perform the chi squared method.
    */
    void CSM();

    /**
     *  @brief  Does the point defined by pfoHCalToHadEnergy and pfoECalToHadEnergy fall within 3 sigma of the ideal distribution
     * 
     *  @param pfoECalToHadEnergy  : Hadronic energy in ECal
     *  @param pfoHCalToHadEnergy  : Hadronic energy in HCal
    */
    bool ThreeSigmaCut(float pfoECalToHadEnergy, float pfoHCalToHadEnergy);

    /**
     *  @brief  Plot m_histogram.
    */
    void Plot();

// Non trivial setting on initialisation
    float               m_hCalToHadResolutionConstant;  ///< HCal coefficient for stoichastic (sigma_E/E ~ coefficeint /sqrt(E)) process default 0.55
    float               m_eCalToHadResolutionConstant;  ///< ECal coefficient for stoichastic (sigma_E/E ~ coefficeint /sqrt(E)) process default 0.55

// Trivial setting on initialisation
    TH2F               *m_histogram;                    ///< 2D histogram of ECalToHad vs HCalToHad energy

typedef std::vector<float> FloatVector;
    FloatVector         m_eCalEnergyAccepted;           ///< Vector of ECalToHad energies passing three sigma cuts
    FloatVector         m_hCalEnergyAccepted;           ///< Vector of HCalToHad energies passing three sigma cuts
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 * 
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  eCalDigitisation to receive the application parameters
 * 
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], ChiSquaredMethod &chiSquaredMethod);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    try
    {
        ChiSquaredMethod chiSquaredMethod;

        if (!ParseCommandLine(argc, argv, chiSquaredMethod))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        chiSquaredMethod.Process();

        std::string dataFileName(chiSquaredMethod.m_outputPath + "Calibration.txt");
        ofstream data_file(dataFileName.c_str(), std::ios_base::app);
        data_file << "_____________________________________________________________________________________" << std::endl;
        data_file << "Hadronic Energy Scale PandoraPFA Calibration performed using the Chi Squared Method." << std::endl << std::endl;
        data_file << "For kaon energy                                    : " << chiSquaredMethod.m_trueEnergy << " : " <<std::endl;
        data_file << "m_eCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_eCalToHadInterceptMinChi2 << " : " <<std::endl;
        data_file << "m_hCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_hCalToHadInterceptMinChi2 << " : " <<std::endl;
        data_file << "m_minChi2                                          : " << chiSquaredMethod.m_minChi2 << " : " <<std::endl;
        data_file << "The total number of events considered was          : " << chiSquaredMethod.m_nEntriesHist << " : " <<std::endl<<std::endl;

        std::cout << "_____________________________________________________________________________________" << std::endl;
        std::cout << "For kaon energy                                    : " << chiSquaredMethod.m_trueEnergy << " : " <<std::endl;
        std::cout << "m_eCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_eCalToHadInterceptMinChi2 << " : " <<std::endl;
        std::cout << "m_hCalToHadInterceptMinChi2                        : " << chiSquaredMethod.m_hCalToHadInterceptMinChi2 << " : " <<std::endl;
        std::cout << "m_minChi2                                          : " << chiSquaredMethod.m_minChi2 << " : " <<std::endl;
        std::cout << "The total number of events considered was          : " << chiSquaredMethod.m_nEntriesHist << " : " <<std::endl<< std::endl;

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

ChiSquaredMethod::ChiSquaredMethod() :
    m_inputKaonLRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_calibrationAccuracy(0.005),
    m_outputPath(""),
    m_numberHCalLayers(48),
    m_eCalToHadInterceptMinChi2(std::numeric_limits<float>::max()),
    m_hCalToHadInterceptMinChi2(std::numeric_limits<float>::max()),
    m_minChi2(std::numeric_limits<float>::max()),
    m_nEntriesHist(0),
    m_hCalToHadResolutionConstant(0.55), // Resolution constants from single particle calibration
    m_eCalToHadResolutionConstant(0.55), // Resolution constants from single particle calibration
    m_histogram(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ChiSquaredMethod::~ChiSquaredMethod()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::Process()
{
    this->CreateHistogram();
    this->FillHistogram();
    this->CSM();
    this->Plot();
}
//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::CreateHistogram()
{
    std::string Name = "histCSM";
    std::string Title = "Hadronic Energy Scale PandoraPFA Calibration, Chi Squared Method";
    int binNumber = static_cast<int>( 1.5 / m_calibrationAccuracy );
    float maxHistogramEnergy = 1.5 * m_trueEnergy;
    m_histogram = new TH2F(Name.c_str(),Title.c_str(),binNumber,0,maxHistogramEnergy,binNumber,0,maxHistogramEnergy);
    m_histogram->GetYaxis()->SetTitle("Hadronic Energy Measured in the ECal / GeV");
    m_histogram->GetXaxis()->SetTitle("Hadronic Energy Measured in the HCal / GeV");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::FillHistogram()
{
    int nPfoTargetsTotal, nPfoTargetsNeutralHadrons, nPfosTotal, nPfosNeutralHadrons, pfoMinHCalLayerToEdge;
    float pfoECalToHadEnergy, pfoHCalToHadEnergy;

    TChain *pTChain = new TChain("PfoAnalysisTree");
    pTChain->Add(m_inputKaonLRootFiles.c_str());

    unsigned int nEntriesHistogram(pTChain->GetEntries());

    pTChain->SetBranchAddress("pfoECalToHadEnergy",&pfoECalToHadEnergy);
    pTChain->SetBranchAddress("pfoHCalToHadEnergy",&pfoHCalToHadEnergy);
    pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    pTChain->SetBranchAddress("nPfoTargetsNeutralHadrons",&nPfoTargetsNeutralHadrons);
    pTChain->SetBranchAddress("nPfosTotal",&nPfosTotal);
    pTChain->SetBranchAddress("nPfosNeutralHadrons",&nPfosNeutralHadrons);
    pTChain->SetBranchAddress("pfoMinHCalLayerToEdge",&pfoMinHCalLayerToEdge);

    int containedLayerExclusion = ceil(m_numberHCalLayers * 0.1);

    for (unsigned int i = 0; i < nEntriesHistogram ; i++) 
    {
        pTChain->GetEntry(i);

        if (ThreeSigmaCut(pfoECalToHadEnergy,pfoHCalToHadEnergy) && nPfoTargetsTotal == 1 && nPfoTargetsNeutralHadrons == 1 && nPfosTotal == 1 && nPfosNeutralHadrons == 1 && 
            pfoMinHCalLayerToEdge > containedLayerExclusion)
        {
            m_eCalEnergyAccepted.push_back(pfoECalToHadEnergy);
            m_hCalEnergyAccepted.push_back(pfoHCalToHadEnergy);
            m_nEntriesHist++;
            m_histogram->Fill(pfoHCalToHadEnergy,pfoECalToHadEnergy);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::CSM()
{
    float perpDist = std::numeric_limits<float>::max();
    float perpDistSigma = std::numeric_limits<float>::max();
    float eCalToHadIntercept = std::numeric_limits<float>::max();
    float hCalToHadIntercept = std::numeric_limits<float>::max();
    float Chi2 = std::numeric_limits<float>::max();

    unsigned int steps = static_cast<int>( 1.f / m_calibrationAccuracy );
    float stepSize = m_trueEnergy * m_calibrationAccuracy;
    float kineticEnergy = m_trueEnergy - 0.497672;

    for (unsigned int p = 0; p < steps; p++)
    {
        eCalToHadIntercept = (kineticEnergy * 0.5) + (stepSize * p );

        for (unsigned int q = 0; q < steps ; q++)
        {
            Chi2 = 0;

            hCalToHadIntercept = (kineticEnergy * 0.5) + (stepSize * q );

            for (unsigned int j = 0; j < m_eCalEnergyAccepted.size() ; j++)
            {
                float eCalSigmaAccepted = m_eCalToHadResolutionConstant * sqrt(m_eCalEnergyAccepted[j]);
                float hCalSigmaAccepted = m_hCalToHadResolutionConstant * sqrt(m_hCalEnergyAccepted[j]);
                
                perpDist = abs( ( (m_hCalEnergyAccepted[j] * eCalToHadIntercept) + (m_eCalEnergyAccepted[j] * hCalToHadIntercept) - (hCalToHadIntercept * eCalToHadIntercept) ) / sqrt( (hCalToHadIntercept * hCalToHadIntercept) + (eCalToHadIntercept * eCalToHadIntercept) ) );
                perpDistSigma = sqrt( ( pow (hCalSigmaAccepted * eCalToHadIntercept, 2) + pow(eCalSigmaAccepted * hCalToHadIntercept, 2) ) / ( pow(hCalToHadIntercept,2) + pow(eCalToHadIntercept,2) ) );

                Chi2 = Chi2 + pow(perpDist/perpDistSigma ,2);
            }

            if (Chi2 < m_minChi2)
            {
                m_eCalToHadInterceptMinChi2 = eCalToHadIntercept;
                m_hCalToHadInterceptMinChi2 = hCalToHadIntercept;
                m_minChi2 = Chi2;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ChiSquaredMethod::ThreeSigmaCut(float pfoECalToHadEnergy, float pfoHCalToHadEnergy)
{
    float idealHCalToHadIntercept = m_trueEnergy - 0.497672; // Want to reconstruct kinetic energy
    float idealECalToHadIntercept = m_trueEnergy - 0.497672; // Want to reconstruct kinetic energy

    float pfoECalToHadEnergySigma = m_eCalToHadResolutionConstant * sqrt(pfoECalToHadEnergy);
    float pfoHCalToHadEnergySigma = m_hCalToHadResolutionConstant * sqrt(pfoHCalToHadEnergy);

    // Formula for perpendicular distance of point m_pfoE/HCalToHad to ideal distribution
    float perpDist = abs( ( (pfoHCalToHadEnergy * idealECalToHadIntercept) + (pfoECalToHadEnergy * idealHCalToHadIntercept) - (idealHCalToHadIntercept * idealECalToHadIntercept) ) / sqrt( (idealHCalToHadIntercept * idealHCalToHadIntercept) + (idealECalToHadIntercept * idealECalToHadIntercept) ) );
    float perpDistSigma = sqrt( ( pow (pfoHCalToHadEnergySigma * idealECalToHadIntercept, 2) + pow(pfoECalToHadEnergySigma * idealHCalToHadIntercept, 2) ) / ( pow(idealHCalToHadIntercept,2) + pow(idealECalToHadIntercept,2) ) );

    if ( perpDist <= 3 * perpDistSigma)
        return true;

    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ChiSquaredMethod::Plot()
{
    TCanvas *pCanvas = new TCanvas("2D Historgram","2D Histogram",200,10,1400,1400);
    pCanvas->cd();
    m_histogram->Draw("COLZ");

    pCanvas->SaveAs(m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Chi_Sqaured_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.png");
    pCanvas->SaveAs(m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Chi_Sqaured_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.C");
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ChiSquaredMethod &chiSquaredMethod)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d:e:f")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            chiSquaredMethod.m_inputKaonLRootFiles = optarg;
            break;
        case 'b':
            chiSquaredMethod.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            chiSquaredMethod.m_calibrationAccuracy = atof(optarg);
            break;
        case 'd':
            chiSquaredMethod.m_outputPath = optarg;
            break;
        case 'e':
            chiSquaredMethod.m_numberHCalLayers = atoi(optarg);
            break;
        case 'f':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)   " << std::endl
                      << "    -b value  (mandatory, true energy of KaonL being used for calibration)                    " << std::endl
                      << "    -c value  (optional, fractional accuracy to calibrate H/ECalToHad to, default 0.005)      " << std::endl
                      << "    -d        (mandatory, output path to send results to)                                     " << std::endl
                      << "    -e value  (mandatory, number of HCal layers in simulation, default 48)                    " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
