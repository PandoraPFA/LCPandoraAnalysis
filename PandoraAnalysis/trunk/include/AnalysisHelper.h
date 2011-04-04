/**
 *  @file   PandoraAnalysis/include/AnalysisHelper.h
 * 
 *  @brief  Header file for the analysis helper class.
 * 
 *  $Log: $
 */

#ifndef ANALYSIS_HELPER_H
#define ANALYSIS_HELPER_H 1

class TH1F;

//------------------------------------------------------------------------------------------------------------------------------------------

namespace pandora_analysis
{

/**
 *  @brief  analysis helper class
 */
class AnalysisHelper
{
public:
    /**
     *  @brief  Calculate performance figures for an energy spectrum provided in form of a root th1f
     * 
     *  @param  pTH1F
     *  @param  sigma
     *  @param  sigmasigma
     *  @param  fixDistributionCentre
     */
    static void CalculatePerformance(const TH1F *const pTH1F, float &sigma, float &sigmasigma, bool fixDistributionCentre = true);
};

} // namespace pandora_analysis

#endif // #ifndef ANALYSIS_HELPER_H
