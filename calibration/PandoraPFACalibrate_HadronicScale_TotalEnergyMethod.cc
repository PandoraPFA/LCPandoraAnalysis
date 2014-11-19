/**
 *  @file   PandoraAnalysis/calibration/PandoraPFACalibrate_HadronicScale_TotalEnergyMethod.cc
 * 
 *  @brief  Calculate E/HCalToHadGeVCalibration using the total energy method on kaonL events
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
#include <fstream>
#include <limits>
#include <stdexcept>

/**
 *  @brief  TotalEnergyMethod class
 */
class TotalEnergyMethod
{
public:
    /**
     *  @brief  Constructor
    */
    TotalEnergyMethod();

    /**
     *  @brief  Destructor
     */
    ~TotalEnergyMethod();

    /**
     *  @brief  Perform the total energy method for hadronic energy scale calibration in PandoraPFA
    */
    void Process();

// Outputs
    float               m_minRMS;                               ///< Narrowest standard deviation found by rescaling E/HCalToHad independantly
    float               m_minRMSECalMultiplier;                 ///< Factor by which ECalToHad energy needs to be scaled by to produce narrowest standard deviation results
    float               m_minRMSHCalMultiplier;                 ///< Factor by which HCalToHad energy needs to be scaled by to produce narrowest standard deviation results
    int                 m_minRMSEvents;                         ///< Number of events in histogram with minimum RMS (different numbers pass 3sigma cut, so not constant)

// Non trivial setting on initialisation
    float               m_trueEnergy;                           ///< Total energy of calibration particle
    std::string         m_inputKaonLRootFiles;                  ///< Input root files for E/HCalToHad Calibration
    std::string         m_outputPath;                           ///< Output path to send plots to
    int                 m_numberHCalLayers;                     ///< Number of layers in the HCal

private:

    /**
     *  @brief  Initial coarse search for min RMS E/HCal multipliers to find min RMS (itter = 0.025)
    */
    void CoarseSearch();

    /**
     *  @brief  Final fine search for min RMS E/HCal multipliers to find min RMS (itter = 0.01)
    */
    void FineSearch();

    /**
     *  @brief  Scan E/HCal multipliers in search of minimum RMS (2D cross param scan, 5 points cross chape step size m_itter)
    */
    void Itterator();

    /**
     *  @brief  Does the point defined by m_pfoHCalToHadEnergy and m_pfoECalToHadEnergy fall within 3 sigme of the ideal distribution
    */
    bool ThreeSigmaCut();

    /**
     *  @brief  Sets m_fitRangeLow, m_fitRangeHigh and m_rMSFitRange based on m_fitPercentage using m_histogram
    */
    void RMSFitPercentageRange();

    /**
     *  @brief  Performs a Gaussian fit to m_histogram and sets m_amplitude, m_mean, m_stdDev and m_chi2.  If m_plot on also plots data.
    */
    void Fit();

// Non trivial setting on initialisation
    float               m_kineticEnergy;                        ///< Kinetic energy of KaonL
    float               m_fitPercentage;                        ///< Percentage of data to fit Gaussian to. Percentage with narrowest range fitted
    float               m_eCalToHadResolutionConstant;          ///< Single particel resolution stoichastic coefficient 0.55 default
    float               m_hCalToHadResolutionConstant;          ///< Single particel resolution stoichastic coefficient 0.55 default

// Trivial setting on initialisation
    float               m_pfoECalToHadEnergy;                   ///< Tree variable.  Hadronic energy in ECal.
    float               m_pfoHCalToHadEnergy;                   ///< Tree variable.  Hadronic energy in HCal.
    float               m_fineGranularityECalMultiplierGuess;   ///< Initial guess for ECal multiplier in FineSearch()
    float               m_fineGranularityHCalMultiplierGuess;   ///< Initial guess for HCal multiplier in FineSearch()

    float               m_itter;                                ///< Step size for scan E/HCalMultiplierGuess space.  0.025 in CoarseSearch() 0.01 in FineSearch()
    float               m_eCalMultiplierRangeMin;               ///< Lower edge of range of ECal multipliers to scan over
    float               m_hCalMultiplierRangeMin;               ///< Lower edge of range of HCal multipliers to scan over
    float               m_eCalMultiplierRangeMax;               ///< Upper edge of range of ECal multipliers to scan over
    float               m_hCalMultiplierRangeMax;               ///< Upper edge of range of HCal multipliers to scan over

    TChain             *m_pTChain;                              ///< Chain of root files
    TH1F               *m_histogram;                            ///< Histogram
    float               m_fitRangeLow;                          ///< Low fit edge of m_fitPercentage of data with min RMS
    float               m_fitRangeHigh;                         ///< High fit edge of m_fitPercentage of data with min RMS
    float               m_rMSFitRange;                          ///< RMS of the m_fitPercentage data
    int                 m_binNumber;                            ///< Number of bins on m_histogram
    float               m_maxHistogramEnergy;                   ///< Max energy on m_histogram
    bool                m_plot;                                 ///< Plot results of fit

    float               m_amplitude;                            ///< Amplitude of Gaussian fit
    float               m_mean;                                 ///< Mean of Gaussian fit
    float               m_stdDev;                               ///< Standard deviation of Gaussian fit
    float               m_chi2;                                 ///< Chi squared of Gaussian fit
    int                 m_nEventsHist;                          ///< Number of events in m_histogram
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
bool ParseCommandLine(int argc, char *argv[], TotalEnergyMethod &totalEnergyMethod);

//------------------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    TApplication *pTApplication = NULL;

    gROOT->SetBatch();

    try
    {
        TotalEnergyMethod totalEnergyMethod;

        if (!ParseCommandLine(argc, argv, totalEnergyMethod))
            return 1;

        pTApplication = new TApplication("MyTest", &argc, argv);

        totalEnergyMethod.Process();

        std::string dataFileName(totalEnergyMethod.m_outputPath + "Calibration.txt");
        ofstream data_file(dataFileName.c_str(), std::ios_base::app);
        data_file << "_____________________________________________________________________________________" << std::endl;
        data_file << "Hadronic Energy Scale PandoraPFA Calibration performed using the Total Energy Method." << std::endl << std::endl;
        data_file << "For kaon energy                                    : " << totalEnergyMethod.m_trueEnergy << " : " <<std::endl;
        data_file << "Minimum_RMS                                        : " << totalEnergyMethod.m_minRMS << " : " <<std::endl;
        data_file << "Minimum_RMS_ECal_Multiplier                        : " << totalEnergyMethod.m_minRMSECalMultiplier << " : " <<std::endl;
        data_file << "Minimum_RMS_HCal_Multiplier                        : " << totalEnergyMethod.m_minRMSHCalMultiplier << " : " <<std::endl;
        data_file << "The total number of events considered was          : " << totalEnergyMethod.m_minRMSEvents << " : " <<std::endl<<std::endl;

        std::cout << "_____________________________________________________________________________________" << std::endl;
        std::cout << "Hadronic Energy Scale PandoraPFA Calibration performed using the Total Energy Method." << std::endl << std::endl;
        std::cout << "For kaon energy                                    : " << totalEnergyMethod.m_trueEnergy << " : " <<std::endl;
        std::cout << "Minimum_RMS                                        : " << totalEnergyMethod.m_minRMS << " : " <<std::endl;
        std::cout << "Minimum_RMS_ECal_Multiplier                        : " << totalEnergyMethod.m_minRMSECalMultiplier << " : " <<std::endl;
        std::cout << "Minimum_RMS_HCal_Multiplier                        : " << totalEnergyMethod.m_minRMSHCalMultiplier << " : " <<std::endl;
        std::cout << "The total number of events considered was          : " << totalEnergyMethod.m_minRMSEvents << " : " <<std::endl<<std::endl;

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

TotalEnergyMethod::TotalEnergyMethod() :
    m_minRMS(std::numeric_limits<float>::max()),
    m_minRMSECalMultiplier(std::numeric_limits<float>::max()),
    m_minRMSHCalMultiplier(std::numeric_limits<float>::max()),
    m_minRMSEvents(std::numeric_limits<int>::max()),
    m_trueEnergy(std::numeric_limits<float>::max()),
    m_inputKaonLRootFiles(""),
    m_outputPath(""),
    m_numberHCalLayers(48),
    m_kineticEnergy(std::numeric_limits<float>::max()),
    m_fitPercentage(90.f),
    m_eCalToHadResolutionConstant(0.55),
    m_hCalToHadResolutionConstant(0.55),
    m_pfoECalToHadEnergy(std::numeric_limits<float>::max()),
    m_pfoHCalToHadEnergy(std::numeric_limits<float>::max()),
    m_fineGranularityECalMultiplierGuess(std::numeric_limits<float>::max()),
    m_fineGranularityHCalMultiplierGuess(std::numeric_limits<float>::max()),
    m_itter(std::numeric_limits<float>::max()),
    m_eCalMultiplierRangeMin(std::numeric_limits<float>::max()),
    m_hCalMultiplierRangeMin(std::numeric_limits<float>::max()),
    m_eCalMultiplierRangeMax(std::numeric_limits<float>::max()),
    m_hCalMultiplierRangeMax(std::numeric_limits<float>::max()),
    m_pTChain(NULL),
    m_histogram(NULL),
    m_fitRangeLow(std::numeric_limits<float>::max()),
    m_fitRangeHigh(std::numeric_limits<float>::max()),
    m_rMSFitRange(std::numeric_limits<float>::max()),
    m_binNumber(100),
    m_maxHistogramEnergy(std::numeric_limits<float>::max()),
    m_plot(false),
    m_amplitude(std::numeric_limits<float>::max()),
    m_mean(std::numeric_limits<float>::max()),
    m_stdDev(std::numeric_limits<float>::max()),
    m_chi2(std::numeric_limits<float>::max()),
    m_nEventsHist(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TotalEnergyMethod::~TotalEnergyMethod()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TotalEnergyMethod::Process()
{
    m_maxHistogramEnergy = 2.5 * m_trueEnergy;
    m_kineticEnergy = m_trueEnergy - 0.497672;
    this->CoarseSearch();
    this->FineSearch();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TotalEnergyMethod::CoarseSearch()
{
    m_itter = 0.025;

    float eCalMultiplierGuess = 1.0;
    float hCalMultiplierGuess = 1.0;

    m_eCalMultiplierRangeMin = eCalMultiplierGuess - m_itter;
    m_hCalMultiplierRangeMin = hCalMultiplierGuess - m_itter;
    m_eCalMultiplierRangeMax = eCalMultiplierGuess + m_itter;
    m_hCalMultiplierRangeMax = hCalMultiplierGuess + m_itter;

    Itterator();

    // 0 is false 1 is true
    Bool_t eCalGuessBig = std::fabs(m_minRMSECalMultiplier - (eCalMultiplierGuess - m_itter)) < m_itter/2.0;
    Bool_t eCalGuessSmall = std::fabs(m_minRMSECalMultiplier - (eCalMultiplierGuess + m_itter)) < m_itter/2.0;
    Bool_t hCalGuessBig = std::fabs(m_minRMSHCalMultiplier - (hCalMultiplierGuess - m_itter)) < m_itter/2.0;
    Bool_t hCalGuessSmall = std::fabs(m_minRMSHCalMultiplier - (hCalMultiplierGuess + m_itter)) < m_itter/2.0;
    
    while (eCalGuessBig == 1 || eCalGuessSmall == 1 || hCalGuessBig == 1 || hCalGuessSmall == 1)
    {
        if (eCalGuessBig == 1)
        {
            eCalMultiplierGuess = eCalMultiplierGuess - m_itter;
        }

        else if (eCalGuessSmall == 1)
        {
            eCalMultiplierGuess = eCalMultiplierGuess + m_itter;
        }

        else if (hCalGuessBig == 1)
        {
            hCalMultiplierGuess = hCalMultiplierGuess - m_itter;
        }

        else if (hCalGuessSmall == 1)
        {
            hCalMultiplierGuess = hCalMultiplierGuess + m_itter;
        }

        m_eCalMultiplierRangeMin = eCalMultiplierGuess - m_itter;
        m_hCalMultiplierRangeMin = hCalMultiplierGuess - m_itter;
        m_eCalMultiplierRangeMax = eCalMultiplierGuess + m_itter;
        m_hCalMultiplierRangeMax = hCalMultiplierGuess + m_itter;

        Itterator();

        eCalGuessBig = std::fabs(m_minRMSECalMultiplier - (eCalMultiplierGuess - m_itter)) < m_itter/2.0;
        eCalGuessSmall = std::fabs(m_minRMSECalMultiplier - (eCalMultiplierGuess + m_itter)) < m_itter/2.0;
        hCalGuessBig = std::fabs(m_minRMSHCalMultiplier - (hCalMultiplierGuess - m_itter)) < m_itter/2.0;
        hCalGuessSmall = std::fabs(m_minRMSHCalMultiplier - (hCalMultiplierGuess + m_itter)) < m_itter/2.0;
    }
    m_fineGranularityECalMultiplierGuess = m_minRMSECalMultiplier;
    m_fineGranularityHCalMultiplierGuess = m_minRMSHCalMultiplier;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TotalEnergyMethod::FineSearch()
{
    m_itter = 0.01;

    float fineGranularityECalMultiplierGuess = m_fineGranularityECalMultiplierGuess;
    float fineGranularityHCalMultiplierGuess = m_fineGranularityHCalMultiplierGuess;

    m_eCalMultiplierRangeMin = fineGranularityECalMultiplierGuess - m_itter;
    m_hCalMultiplierRangeMin = fineGranularityHCalMultiplierGuess - m_itter;
    m_eCalMultiplierRangeMax = fineGranularityECalMultiplierGuess + m_itter;
    m_hCalMultiplierRangeMax = fineGranularityHCalMultiplierGuess + m_itter;

    Itterator();

    Bool_t eCalGuessBig = std::fabs(m_minRMSECalMultiplier - (fineGranularityECalMultiplierGuess - m_itter)) < m_itter/2.0;
    Bool_t eCalGuessSmall = std::fabs(m_minRMSECalMultiplier - (fineGranularityECalMultiplierGuess + m_itter)) < m_itter/2.0;
    Bool_t hCalGuessBig = std::fabs(m_minRMSHCalMultiplier - (fineGranularityHCalMultiplierGuess - m_itter)) < m_itter/2.0;
    Bool_t hCalGuessSmall = std::fabs(m_minRMSHCalMultiplier - (fineGranularityHCalMultiplierGuess + m_itter)) < m_itter/2.0;

    while (eCalGuessBig == 1 || eCalGuessSmall == 1 || hCalGuessBig == 1 || hCalGuessSmall == 1)
    {
        if (eCalGuessBig == 1)
        {
            fineGranularityECalMultiplierGuess = fineGranularityECalMultiplierGuess - m_itter;
        }

        else if (eCalGuessSmall == 1)
        {
            fineGranularityECalMultiplierGuess = fineGranularityECalMultiplierGuess + m_itter;
        }

        else if (hCalGuessBig == 1)
        {
            fineGranularityHCalMultiplierGuess = fineGranularityHCalMultiplierGuess - m_itter;
        }

        else if (hCalGuessSmall == 1)
        {
            fineGranularityHCalMultiplierGuess = fineGranularityHCalMultiplierGuess + m_itter;
        }

        m_eCalMultiplierRangeMin = fineGranularityECalMultiplierGuess - m_itter;
        m_hCalMultiplierRangeMin = fineGranularityHCalMultiplierGuess - m_itter;
        m_eCalMultiplierRangeMax = fineGranularityECalMultiplierGuess + m_itter;
        m_hCalMultiplierRangeMax = fineGranularityHCalMultiplierGuess + m_itter;

        Itterator();

        eCalGuessBig = std::fabs(m_minRMSECalMultiplier - (fineGranularityECalMultiplierGuess - m_itter)) < m_itter/2.0;
        eCalGuessSmall = std::fabs(m_minRMSECalMultiplier - (fineGranularityECalMultiplierGuess + m_itter)) < m_itter/2.0;
        hCalGuessBig = std::fabs(m_minRMSHCalMultiplier - (fineGranularityHCalMultiplierGuess - m_itter)) < m_itter/2.0;
        hCalGuessSmall = std::fabs(m_minRMSHCalMultiplier - (fineGranularityHCalMultiplierGuess + m_itter)) < m_itter/2.0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TotalEnergyMethod::Itterator()
{
    m_pTChain = NULL;
    m_pTChain = new TChain("PfoAnalysisTree");
    m_pTChain->Add(m_inputKaonLRootFiles.c_str());

    int nPfoTargetsTotal;
    int nPfoTargetsNeutralHadrons;
    int nPfosTotal;
    int nPfosNeutralHadrons;
    int pfoMinHCalLayerToEdge;
    const unsigned int nEntries = m_pTChain->GetEntries();

    m_pTChain->SetBranchAddress("pfoECalToHadEnergy",&m_pfoECalToHadEnergy);
    m_pTChain->SetBranchAddress("pfoHCalToHadEnergy",&m_pfoHCalToHadEnergy);
    m_pTChain->SetBranchAddress("nPfoTargetsTotal",&nPfoTargetsTotal);
    m_pTChain->SetBranchAddress("nPfoTargetsNeutralHadrons",&nPfoTargetsNeutralHadrons);
    m_pTChain->SetBranchAddress("nPfosTotal",&nPfosTotal);
    m_pTChain->SetBranchAddress("nPfosNeutralHadrons",&nPfosNeutralHadrons);
    m_pTChain->SetBranchAddress("pfoMinHCalLayerToEdge",&pfoMinHCalLayerToEdge);

    float eCalMultiplierAverage = (m_eCalMultiplierRangeMin + m_eCalMultiplierRangeMax)/2.0; 

    int containedLayerExclusion = ceil(m_numberHCalLayers * 0.1);

    for (float hCalMultiplier = m_hCalMultiplierRangeMin; hCalMultiplier <= m_hCalMultiplierRangeMax; hCalMultiplier+=m_itter)
    {
        m_histogram = new TH1F("KaonL Total Energy Histogram", "KaonL Total Energy Histogram", m_binNumber, 0, m_maxHistogramEnergy);

        for (unsigned int k = 0; k < nEntries; k++) 
        {
            m_pTChain->GetEntry(k);

            if (ThreeSigmaCut() && 1 == nPfoTargetsTotal && 1 == nPfoTargetsNeutralHadrons && 1 == nPfosTotal && 1 == nPfosNeutralHadrons && containedLayerExclusion < pfoMinHCalLayerToEdge)
                m_histogram->Fill((eCalMultiplierAverage*m_pfoECalToHadEnergy) + (hCalMultiplier*m_pfoHCalToHadEnergy));
        }

        RMSFitPercentageRange();
        Fit();

        if (m_stdDev < m_minRMS)
        {
            m_minRMS = m_stdDev;
            m_minRMSECalMultiplier = eCalMultiplierAverage;
            m_minRMSHCalMultiplier = hCalMultiplier;
            m_minRMSEvents = m_histogram->GetEntries();

            m_plot = true;
            Fit();
            m_plot = false;
        }
        m_histogram = NULL;
    }

    float hCalMultiplierAverage = (m_hCalMultiplierRangeMin + m_hCalMultiplierRangeMax)/2.0; 

// m_itter is 2*m_itter as hCalMultiplierAverage and eCalMultiplierAverage covered in first loop
    for (float eCalMultiplier = m_eCalMultiplierRangeMin; eCalMultiplier <= m_eCalMultiplierRangeMax; eCalMultiplier+=(2*m_itter))
    {
        m_histogram = new TH1F("KaonL Total Energy Histogram", "KaonL Total Energy Histogram", m_binNumber, 0, m_maxHistogramEnergy);
        
        for (unsigned int k = 0; k < nEntries ; k++) 
        {
            m_pTChain->GetEntry(k);

            if (ThreeSigmaCut() && 1 == nPfoTargetsTotal && 1 == nPfoTargetsNeutralHadrons && 1 == nPfosTotal && 1 == nPfosNeutralHadrons && containedLayerExclusion < pfoMinHCalLayerToEdge)
                m_histogram->Fill((eCalMultiplier*m_pfoECalToHadEnergy) + (hCalMultiplierAverage*m_pfoHCalToHadEnergy));
        }

        RMSFitPercentageRange();
        Fit();

        if (m_stdDev < m_minRMS)
        {
            m_minRMS = m_stdDev;
            m_minRMSECalMultiplier = eCalMultiplier;
            m_minRMSHCalMultiplier = hCalMultiplierAverage;
            m_minRMSEvents = m_histogram->GetEntries();

            m_plot = true;
            Fit();
            m_plot = false;
        }
        m_histogram = NULL;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TotalEnergyMethod::ThreeSigmaCut()
{
    float idealHCalToHadIntercept = m_kineticEnergy;
    float idealECalToHadIntercept = m_kineticEnergy;

    float pfoECalToHadEnergySigma = m_eCalToHadResolutionConstant * sqrt(m_pfoECalToHadEnergy);
    float pfoHCalToHadEnergySigma = m_hCalToHadResolutionConstant * sqrt(m_pfoHCalToHadEnergy);

    // Formula for perpendicular distance of point m_pfoE/HCalToHad to ideal distribution
    float perpDist = abs( ( (m_pfoHCalToHadEnergy * idealECalToHadIntercept) + (m_pfoECalToHadEnergy * idealHCalToHadIntercept) - (idealHCalToHadIntercept * idealECalToHadIntercept) ) / sqrt( (idealHCalToHadIntercept * idealHCalToHadIntercept) + (idealECalToHadIntercept * idealECalToHadIntercept) ) );
    float perpDistSigma = sqrt( ( pow (pfoHCalToHadEnergySigma * idealECalToHadIntercept, 2) + pow(pfoECalToHadEnergySigma * idealHCalToHadIntercept, 2) ) / ( pow(idealHCalToHadIntercept,2) + pow(idealECalToHadIntercept,2) ) );

    if ( perpDist <= 3 * perpDistSigma)
        return true;

    else 
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TotalEnergyMethod::RMSFitPercentageRange()
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

void TotalEnergyMethod::Fit()
{
    try
    {
        if (m_histogram->GetEntries() != 0)
        {
            // Fit Function
            TF1 *gaussianFitFunc = new TF1("gaussianFitFunc","[0] * TMath::Exp( -0.5 * [2] * TMath::Power(x-[1],2) )");

            // Initial Param Guess
            float histAmp = m_histogram->GetBinContent(m_histogram->GetMaximumBin());
            float histMean = m_histogram->GetMean();
            float histRMS = std::pow(m_histogram->GetRMS(),-2);

            gaussianFitFunc->SetParameters(0, histAmp);
            gaussianFitFunc->SetParameters(1, histMean);
            gaussianFitFunc->SetParameters(2, histRMS);
            gaussianFitFunc->SetParLimits(2,0,100);

            // Perform Fit and Get Results
            TFitResultPtr pTFitResultPtr = m_histogram->Fit(gaussianFitFunc, "SM", "", m_fitRangeLow, m_fitRangeHigh);

            bool isValidFit = pTFitResultPtr->IsValid();
            int fitQuality = pTFitResultPtr->CovMatrixStatus();
            int count = 0;

            //std::cout << "isValidFit          :" << isValidFit << std::endl;
            //std::cout << "fitQuality          :" << fitQuality << std::endl;
            //std::cout << "Count               :" << count << std::endl;

            // While loop until fit converges
            while(1 != isValidFit || 3 != fitQuality || gaussianFitFunc->GetParameter(1) < m_fitRangeLow || gaussianFitFunc->GetParameter(1) > m_fitRangeHigh)
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

                gaussianFitFunc->SetParameters(0, RescaleFactor1 * histAmp);
                gaussianFitFunc->SetParameters(1, RescaleFactor2 * histMean);
                gaussianFitFunc->SetParameters(2, RescaleFactor3 * histRMS);

                gaussianFitFunc->SetParLimits(0, LowerFraction1 * RescaleFactor1 * histAmp, UpperFraction1 * RescaleFactor1 * histAmp);
                gaussianFitFunc->SetParLimits(1, LowerFraction2 * RescaleFactor2 * histMean, UpperFraction2 * RescaleFactor2 * histMean);
                gaussianFitFunc->SetParLimits(2, LowerFraction3 * RescaleFactor3 * histRMS, UpperFraction3 * RescaleFactor3 * histRMS);

                pTFitResultPtr = m_histogram->Fit(gaussianFitFunc, "SM", "", m_fitRangeLow, m_fitRangeHigh);

                isValidFit = pTFitResultPtr->IsValid();
                fitQuality = pTFitResultPtr->CovMatrixStatus();
                count++;

                //std::cout << "isValidFit          :" << isValidFit << std::endl;
                //std::cout << "fitQuality          :" << fitQuality << std::endl;
                //std::cout << "Count               :" << count << std::endl;

                if (count > 1000)
                    throw std::runtime_error("Fitting attempts exceeding maximum limit (1000).");
            }

            m_amplitude = gaussianFitFunc->GetParameter(0);
            m_mean = gaussianFitFunc->GetParameter(1);
            m_stdDev = std::pow(gaussianFitFunc->GetParameter(2),-0.5);
            m_chi2 = gaussianFitFunc->GetChisquare();

            if (m_plot)
            {
                std::string canvasName = "histTEM";
                std::string canvasTitle = "Hadronic Energy Scale PandoraPFA Calibration, Total Energy Method";
                TCanvas *pCanvas = new TCanvas(canvasName.c_str(),canvasTitle.c_str(),200,10,600,500);
                pCanvas->cd();
                m_histogram->GetYaxis()->SetTitle("Entries");
                m_histogram->GetXaxis()->SetTitle("Scaled Hadronic PFO Energy / GeV");
                m_histogram->Draw("");
                gaussianFitFunc->Draw("same");

                TString pngOutputFilename = m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Total_Energy_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.png";
                TString dotCOutputFilename = m_outputPath + "PandoraPFA_Calibration_Hadronic_Energy_Scale_Total_Energy_Method_" + TString::Format("%i",int(m_trueEnergy + 0.5)) + "_GeV_KaonL.C";

                pCanvas->SaveAs(pngOutputFilename);
                pCanvas->SaveAs(dotCOutputFilename);
                delete pCanvas;
            }
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

bool ParseCommandLine(int argc, char *argv[], TotalEnergyMethod &totalEnergyMethod)
{
    int c(0);

    while (((c = getopt(argc, argv, "e:i:g:o:f")) != -1) || (argc == 1))
    {
        switch (c)
        {
        case 'e':
            totalEnergyMethod.m_trueEnergy = atof(optarg);
            break;
        case 'i':
            totalEnergyMethod.m_inputKaonLRootFiles = optarg;
            break;
        case 'g':
            totalEnergyMethod.m_numberHCalLayers = atoi(optarg);
            break;
        case 'o':
            totalEnergyMethod.m_outputPath = optarg;
            break;
        case 'f':
        default:
            std::cout << std::endl << "Calibrate " << std::endl
                      << "    -e value  (mandatory, true energy of KaonL being used for calibration)" << std::endl
                      << "    -i        (mandatory, input file name(s), can include wildcards if string is in quotes)" << std::endl
                      << "    -g value  (mandatory, number of HCal layers in simulation, default 48)" << std::endl
                      << "    -o value  (mandatory, output path to send results to)" << std::endl
                      << std::endl;
            return false;
        }
    }
    return true;
}