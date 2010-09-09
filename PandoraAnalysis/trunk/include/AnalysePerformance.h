/**
 *  @file   PandoraAnalysis/include/AnalysePerformance.h
 * 
 *  @brief  Header file for pandora analyse performance binary.
 * 
 *  $Log: $
 */

#ifndef ANALYSE_PERFORMANCE_H
#define ANALYSE_PERFORMANCE_H 1

class TH1F;
class TRandom;
class TFile;

/**
 *  @brief  Analyse histograms in specified root file
 * 
 *  @param  pTFile address of the root file
 *  @param  outputRootFileName the output root file name
 */
void AnalyseHistograms(TFile *pTFile, const std::string &outputRootFileName);

/**
 *  @brief  Calculate rms for distribution is specified histogram
 * 
 *  @param  pTH1F address of the root histogram
 *  @param  fixDistributionCentre 
 */
void CalcRms(TH1 *pTH1F, bool fixDistributionCentre = true);

float           m_mean;             ///< 
float           m_rms;              ///< 
float           m_low;              ///< 
float           m_high;             ///< 
float           m_sigma;            ///< 
float           m_sigmasigma;       ///< 

TRandom        *m_pRandom;          ///< 

#endif // #ifndef ANALYSE_PERFORMANCE_H
