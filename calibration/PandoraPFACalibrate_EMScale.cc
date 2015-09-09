/**
 *  @file   PandoraAnalysis/calibration/PandoraPFACalibrate_EMScale.cc
 * 
 *  @brief  Gaussian fit to PFO energy for photon events.  Used for setting EM scale in PandoraPFA 
 *          (ECalToEMGeVCalibration/HCalToEMGeVCalibration).
 * 
 *  $Log: $
 */

#include "TApplication.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TTree.h"

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  ECalToEM class
 */
class ECalToEM 
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pFitPercentage Percentage of data to fit Gaussian to. Percentage with narrowest range fitted
    */
    ECalToEM();

    /**
     *  @brief  Destructor
     */
    ~ECalToEM();

    /**
     *  @brief  Create histogram, find the limits for the percentage fit and perform Gaussian fit over this range
    */
    void Process();

// Non trivial setting on initialisation
    std::string     m_inputPhotonRootFiles; ///< Input root files - Photons
    float           m_trueEnergy;           ///< True energy of photons being simulated
    float           m_calibrationAccuracy;  ///< Fractional accuracy to calibrate ECalToEM to
    std::string     m_outputPath;           ///< Output path to send plots to
    float           m_fitPercentage;        ///< Percentage of data to fit Gaussian to. Percentage with narrowest range fitted

// Outputs
    float           m_amplitude;            ///< Amplitude of Gaussian fit
    float           m_mean;                 ///< Mean of Gaussian fit
    float           m_stdDev;               ///< Standard deviation of Gaussian fit
    int             m_nEventsECalHist;      ///< Number of events in m_histogram

private:
    /**
     *  @brief  Set m_maxHistogramEnergy
    */
    void PrepareHistogram();

    /**
     *  @brief  Create ECal calo hit energy histrogram (m_histogram)
    */
    void CreateHistogram();

    /**
     *  @brief  Fill ECal calo hit energy histrogram 
    */
    void FillHistogram();

    /**
     *  @brief  Set fitRangeLow, fitRangeHigh and RMSFitPercentage
    */
    void RMSFitPercentageRange();

    /**
     *  @brief  Perform Gaussian fit to histogram
    */
    void Fit();

// Trivial setting on initialisation
    TChain *m_pTChain;              ///< Chain of root files
    TH1F   *m_histogram;            ///< Histogram
    float   m_fitRangeLow;          ///< Low fit edge
    float   m_fitRangeHigh;         ///< High fit edge
    float   m_rMSFitRange;          ///< RMS of the fitRange data
    float   m_maxHistogramEnergy;   ///< Max energy on m_histogram
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
bool ParseCommandLine(int argc, char *argv[], ECalToEM &eCalToEM);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        ECalToEM eCalToEM;

        if (!ParseCommandLine(argc, argv, eCalToEM))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        eCalToEM.Process();

        std::string dataFileName(eCalToEM.m_outputPath+"Calibration.txt");
        ofstream data_file(dataFileName.c_str(), std::ios_base::app);

        data_file << "_____________________________________________________________________________________" << std::endl;
        data_file << "Electromagnetic Energy Scale PandoraPFA Calibration." << std::endl << std::endl;
        data_file << "For photon energy                                  : " << eCalToEM.m_trueEnergy << " : " <<std::endl;
        data_file << "Gaussian fit to total PFO energy for single photon : " << std::endl;
        data_file << "reconstructed as single photon has the following   : " << std::endl;
        data_file << "parameters, data used 90% min RMS:                 : " << std::endl;
        data_file << "Amplitude                                          : " << eCalToEM.m_amplitude << " : " <<std::endl;
        data_file << "ECalToEM Mean                                      : " << eCalToEM.m_mean << " : " <<std::endl;
        data_file << "Standard Deviation                                 : " << eCalToEM.m_stdDev << " : " <<std::endl;
        data_file << "The total number of events considered was          : " << eCalToEM.m_nEventsECalHist << " : " <<std::endl<<std::endl;

        std::cout << "_____________________________________________________________________________________" << std::endl;
        std::cout << "Electromagnetic Energy Scale PandoraPFA Calibration." << std::endl << std::endl;
        std::cout << "For photon energy                                  : " << eCalToEM.m_trueEnergy << " : " <<std::endl;
        std::cout << "Gaussian fit to total PFO energy for single photon : " << std::endl;
        std::cout << "reconstructed as single photon has the following   : " << std::endl;
        std::cout << "parameters, data used 90% min RMS:                 : " << std::endl;
        std::cout << "Amplitude                                          : " << eCalToEM.m_amplitude << " : " <<std::endl;
        std::cout << "ECalToEM Mean                                      : " << eCalToEM.m_mean << " : " <<std::endl;
        std::cout << "Standard Deviation                                 : " << eCalToEM.m_stdDev << " : " <<std::endl;
        std::cout << "The total number of events considered was          : " << eCalToEM.m_nEventsECalHist << " : " <<std::endl<<std::endl;

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

ECalToEM::ECalToEM() :
    m_inputPhotonRootFiles(""),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_calibrationAccuracy(0.005),
    m_outputPath(""),
    m_fitPercentage(90.f),
    m_amplitude(std::numeric_limits<float>::max()),
    m_mean(std::numeric_limits<float>::max()),
    m_stdDev(std::numeric_limits<float>::max()),
    m_nEventsECalHist(0),
    m_pTChain(NULL),
    m_histogram(NULL),
    m_fitRangeLow(std::numeric_limits<float>::max()),
    m_fitRangeHigh(std::numeric_limits<float>::max()),
    m_rMSFitRange(std::numeric_limits<float>::max()),
    m_maxHistogramEnergy(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ECalToEM::~ECalToEM()
{
}

//--------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::Process()
{
    m_pTChain = new TChain("PfoAnalysisTree");
    m_pTChain->Add(m_inputPhotonRootFiles.c_str());
    this->PrepareHistogram();
    this->CreateHistogram();
    this->FillHistogram();
    this->RMSFitPercentageRange();
    this->Fit();
}

//--------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::PrepareHistogram()
{
    float pfoECalToEmEnergy, pfoHCalToEmEnergy;
    int nPfoTargetsTotal, nPfoTargetsPhotons, nPfosTotal, nPfosPhotons;

    m_pTChain->SetBranchAddress("pfoECalToEmEnergy",&pfoECalToEmEnergy);
    m_pTChain->SetBranchAddress("pfoHCalToEmEnergy",&pfoHCalToEmEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsPhotons",&nPfoTargetsPhotons);
    m_pTChain->SetBranchAddress("nPfosTotal",&nPfosTotal);
    m_pTChain->SetBranchAddress("nPfosPhotons",&nPfosPhotons);

    for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++) 
    {
        m_pTChain->GetEntry(i);

        bool isPrimaryECal(pfoHCalToEmEnergy < 0.01 * m_trueEnergy);
        if (1 == nPfoTargetsTotal && 1 == nPfoTargetsPhotons && 1 == nPfosTotal && 1 == nPfosPhotons && isPrimaryECal)
        {
            if (pfoECalToEmEnergy > m_maxHistogramEnergy)
                m_maxHistogramEnergy = pfoECalToEmEnergy;

            m_nEventsECalHist++;
        }
    }
}

//--------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::CreateHistogram()
{
    std::string Name = "EMEnergyScalePandoraPFA";
    std::string Title = "Electromagentic Energy Scale PandoraPFA Calibration";
    int binNumber = static_cast<int>( m_maxHistogramEnergy / ( m_calibrationAccuracy * m_trueEnergy ) );
    m_histogram = new TH1F(Name.c_str(), Title.c_str(), binNumber, 0., m_maxHistogramEnergy);
    m_histogram->GetXaxis()->SetTitle("Calorimeter Hit Energy / GeV");
    m_histogram->GetYaxis()->SetTitle("Entries");
}

//--------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::FillHistogram()
{
    float pfoECalToEmEnergy, pfoHCalToEmEnergy;
    int nPfoTargetsTotal, nPfoTargetsPhotons, nPfosTotal, nPfosPhotons;

    m_pTChain->SetBranchAddress("pfoECalToEmEnergy",&pfoECalToEmEnergy);
    m_pTChain->SetBranchAddress("pfoHCalToEmEnergy",&pfoHCalToEmEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsPhotons",&nPfoTargetsPhotons);
    m_pTChain->SetBranchAddress("nPfosTotal",&nPfosTotal);
    m_pTChain->SetBranchAddress("nPfosPhotons",&nPfosPhotons);

    for (unsigned int i = 0; i < m_pTChain->GetEntries(); i++) 
    {
        m_pTChain->GetEntry(i);

        bool isPrimaryECal(pfoHCalToEmEnergy < 0.01 * m_trueEnergy);
        if (1 == nPfoTargetsTotal && 1 == nPfoTargetsPhotons && 1 == nPfosTotal && 1 == nPfosPhotons && isPrimaryECal)
            m_histogram->Fill(pfoECalToEmEnergy);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::RMSFitPercentageRange()
{
    static const float FLOAT_MAX(std::numeric_limits<float>::max());

    if (NULL == m_histogram)
        return;

    if (5 > m_histogram->GetEntries())
    {
        std::cout << m_histogram->GetName() << " (" << m_histogram->GetEntries() << " entries) - skipped" << std::endl;
        return;
    }

    // Calculate raw properties of distribution (ie rms100)
    float sum = 0., total = 0.;
    double sx = 0., sxx = 0.;
    const unsigned int nbins(m_histogram->GetNbinsX());

    for (unsigned int i = 0; i <= nbins; ++i)
    {
        const float binx(m_histogram->GetBinLowEdge(i) + (0.5 * m_histogram->GetBinWidth(i)));
        const float yi(m_histogram->GetBinContent(i));
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
    }

    const float rawMean(sx / total);
    const float rawMeanSquared(sxx / total);
    const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));

    sum = 0.;
    unsigned int is0 = 0;

    //  The /10 comes from the fact that for rms 90 the start point for the fit must occur in the first 10% of the data.
    float frac = (1 - (m_fitPercentage/100.0));
    for (unsigned int i = 0; (i <= nbins) && (sum < total * frac); ++i)
    {
        sum += m_histogram->GetBinContent(i);
        is0 = i;
    }

    // Calculate truncated properties
    float rmsmin(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX);
    float high(0.f);

    for (unsigned int istart = 0; istart <= is0; ++istart)
    {
        double sumn = 0.;
        double csum = 0.;
        double sumx = 0.;
        double sumxx = 0.;
        unsigned int iend = 0;

        for (unsigned int i = istart; (i <= nbins) && (csum < (m_fitPercentage/100) * total); ++i)
        {
            const float binx(m_histogram->GetBinLowEdge(i) + (0.5 * m_histogram->GetBinWidth(i)));
            const float yi(m_histogram->GetBinContent(i));
            //csum is the sum of yi from istart and is used to stop the sum when this exceeds X% of data.
            csum += yi;

            if (sumn < (m_fitPercentage/100) * total)
            {
                // These variables define the final sums required once we have considered X% of data, anything else is 
                // continuously overwritten.
                sumn += yi;
                sumx += yi * binx;
                sumxx+= yi * binx * binx;
                iend = i;
            }
        }

        const float localMean(sumx / sumn);
        const float localMeanSquared(sumxx / sumn);
        // Standard deviation formula
        const float localRms(std::sqrt(localMeanSquared - localMean * localMean));

        if (localRms < rmsmin)
        {
            mean = localMean;
            if (istart==0)
            {
                low = 0;
                m_fitRangeLow = 0;
            }
            else
            {
                low = m_histogram->GetBinLowEdge(istart);
                m_fitRangeLow=m_histogram->GetBinLowEdge(istart) + (0.5 * m_histogram->GetBinWidth(istart));
            }
            
            high = m_histogram->GetBinLowEdge(iend);
            rmsmin = localRms;
            m_fitRangeHigh=m_histogram->GetBinLowEdge(iend) + (0.5 * m_histogram->GetBinWidth(iend));
        }
    }
    
    m_rMSFitRange = rmsmin;
    
    std::cout << m_histogram->GetName() << " (" << m_histogram->GetEntries() << " entries), rawrms: " << rawRms << ", rmsx: " << rmsmin
              << " (" << low << "-" << high << "), low_fit and high_fit " << " (" << m_fitRangeLow << "-" << m_fitRangeHigh 
              << "), << mean: " << mean << std::endl;
}

//--------------------------------------------------------------------------------------------------------------------------------------

void ECalToEM::Fit()
{
    TH1F *pTH1F = const_cast<TH1F *>(m_histogram);
    try
    {
        if (m_histogram->GetEntries() != 0)
        {
            // Fit Function
            TF1 *Gaussian_Fit_Func = new TF1("Gaussian_Fit_Func","[0] * TMath::Exp( -0.5 * [2] * TMath::Power(x-[1],2) )");

            // Initial Param Guess
            float Hist_Amp = pTH1F->GetBinContent(pTH1F->GetMaximumBin());
            float Hist_Mean = pTH1F->GetMean();
            float Hist_RMS = std::pow(pTH1F->GetRMS(),-2);

            Gaussian_Fit_Func->SetParameters(0, Hist_Amp);
            Gaussian_Fit_Func->SetParameters(1, Hist_Mean);
            Gaussian_Fit_Func->SetParameters(2, Hist_RMS);
            Gaussian_Fit_Func->SetParLimits(2,0,100);

            // Perform Fit and Get Results
            TFitResultPtr pTFitResultPtr = pTH1F->Fit(Gaussian_Fit_Func, "SM", "", m_fitRangeLow, m_fitRangeHigh);

            bool IsValidFit = pTFitResultPtr->IsValid();
            int FitQuality = pTFitResultPtr->CovMatrixStatus();
            int count = 0;

            // While loop until fit converges
            while (
            IsValidFit != 1 ||
            FitQuality != 3 ||
            Gaussian_Fit_Func->GetParameter(1) < m_fitRangeLow ||
            Gaussian_Fit_Func->GetParameter(1) > m_fitRangeHigh
            )
            {
                const float sign1(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);
                const float sign2(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);
                const float sign3(((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) > 0.5f) ? 1.f : -1.f);

                float RescaleFactor1 = 1 + (0.5 * sign1 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float RescaleFactor2 = 1 + (0.5 * sign2 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float RescaleFactor3 = 1 + (0.5 * sign3 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                float LowerFraction1 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float LowerFraction2 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float LowerFraction3 = 1 - std::numeric_limits<float>::epsilon() - (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                float UpperFraction1 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float UpperFraction2 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));
                float UpperFraction3 = 1 + std::numeric_limits<float>::epsilon() + (0.5 * static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX));

                Gaussian_Fit_Func->SetParameters(0, RescaleFactor1 * Hist_Amp);
                Gaussian_Fit_Func->SetParameters(1, RescaleFactor2 * Hist_Mean);
                Gaussian_Fit_Func->SetParameters(2, RescaleFactor3 * Hist_RMS);

                Gaussian_Fit_Func->SetParLimits(0, LowerFraction1 * RescaleFactor1 * Hist_Amp, UpperFraction1 * RescaleFactor1 * Hist_Amp);
                Gaussian_Fit_Func->SetParLimits(1, LowerFraction2 * RescaleFactor2 * Hist_Mean, UpperFraction2 * RescaleFactor2 * Hist_Mean);
                Gaussian_Fit_Func->SetParLimits(2, LowerFraction3 * RescaleFactor3 * Hist_RMS, UpperFraction3 * RescaleFactor3 * Hist_RMS);

                pTFitResultPtr = pTH1F->Fit(Gaussian_Fit_Func, "SM", "", m_fitRangeLow, m_fitRangeHigh);

                IsValidFit = pTFitResultPtr->IsValid();
                FitQuality = pTFitResultPtr->CovMatrixStatus();
                count++;

                std::cout << "IsValidFit          :" << IsValidFit << std::endl;
                std::cout << "FitQuality          :" << FitQuality << std::endl;
                std::cout << "Count               :" << count << std::endl;

                if (count > 1000)
                    throw std::runtime_error("Fitting attempts exceeding maximum limit (1000).");
            }

            m_amplitude = Gaussian_Fit_Func->GetParameter(0);
            m_mean = Gaussian_Fit_Func->GetParameter(1);
            m_stdDev = std::pow(Gaussian_Fit_Func->GetParameter(2),-0.5);

            std::string canvasName = "ECalDigitisation";
            std::string canvasTitle = "ECal Digitisation";
            TCanvas *pCanvas = new TCanvas(canvasName.c_str(),canvasTitle.c_str(),200,10,600,500);
            pCanvas->cd();
            pTH1F->GetYaxis()->SetTitle("Entries");
            std::string xAxisTitle = "Calorimeter Hit Energy in ECal / GeV";
            pTH1F->GetXaxis()->SetTitle(xAxisTitle.c_str());
            pTH1F->Draw("");
            Gaussian_Fit_Func->Draw("same");

            TString pngOutputFilename = m_outputPath + "PandoraPFA_Calibration_Electromagnetic_Energy_Scale_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Photons.png";
            TString dotCOutputFilename = m_outputPath + "PandoraPFA_Calibration_Electromagnetic_Energy_Scale_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_Photons.C";

            pCanvas->SaveAs(pngOutputFilename);
            pCanvas->SaveAs(dotCOutputFilename);
            delete pCanvas;
        }
        
        else
        {
            std::cout << "Histogram is empty, no fit possible." << std::endl;
        }
    }
    
    catch (const std::runtime_error& e)
    {
        std::cout << "An exception occurred. " << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], ECalToEM &eCalToEM)
{
    int c(0);

    while (((c = getopt(argc, argv, "a:b:c:d:e:f")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'a':
            eCalToEM.m_inputPhotonRootFiles = optarg;
            break;
        case 'b':
            eCalToEM.m_trueEnergy = atof(optarg);
            break;
        case 'c':
            eCalToEM.m_calibrationAccuracy = atof(optarg);
            break;
        case 'd':
            eCalToEM.m_outputPath = optarg;
            break;
        case 'e':
            eCalToEM.m_fitPercentage = atof(optarg);
            break;
        case 'f':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -a        (mandatory, input file name(s), can include wildcards if string is in quotes)           " << std::endl
                      << "    -b value  (mandatory, true energy of photons being used for calibration)                          " << std::endl
                      << "    -c value  (optional, fractional accuracy to calibrate CalibrECAL to, default 0.005)               " << std::endl
                      << "    -d value  (mandatory, output path to send results to)                                             " << std::endl
                      << "    -e value  (optional, fit percentage used for calibration, default 90% of data with narrowest rms) " << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}
