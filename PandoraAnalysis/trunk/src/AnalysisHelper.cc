/**
 *  @file   PandoraAnalysis/src/AnalysisHelper.cc
 * 
 *  @brief  Implementation of the analysis helper class.
 * 
 *  $Log: $
 */

#include "TH1F.h"

#include "AnalysisHelper.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace pandora_analysis
{

/**
 *  @brief  InvalidEnergyException class
 */
class InvalidEnergyException : public std::exception
{
public:
    const char *what() const throw();
};

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisHelper::CalculatePerformance(const TH1F *const pTH1F, float &resolution, float &resolutionError, bool fixDistributionCentre, bool print)
{
    static const float FLOAT_MAX(std::numeric_limits<float>::max());

    if (NULL == pTH1F)
        return;

    if (5 > pTH1F->GetEntries())
    {
        std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
        return;
    }

    // Calculate raw properties of distribution
    float sum = 0., total = 0.;
    double sx = 0., sxx = 0.;
    const unsigned int nbins(pTH1F->GetNbinsX());

    for (unsigned int i = 0; i <= nbins; ++i)
    {
        const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
        const float yi(pTH1F->GetBinContent(i));
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
    }

    const float rawMean(sx / total);
    const float rawMeanSquared(sxx / total);
    const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));

    sum = 0.;
    unsigned int is0 = 0;

    for (unsigned int i = 0; (i <= nbins) && (sum < total / 10.); ++i)
    {
        sum += pTH1F->GetBinContent(i);
        is0 = i;
    }

    // Calculate truncated properties
    float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX), rms(FLOAT_MAX);
    float high(0.f);

    for (unsigned int istart = 0; istart <= is0; ++istart)
    {
        double sumn = 0.;
        double csum = 0.;
        double sumx = 0.;
        double sumxx = 0.;
        unsigned int iend = 0;

        for (unsigned int i = istart; (i <= nbins) && (csum < 0.9 * total); ++i)
        {
            const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
            const float yi(pTH1F->GetBinContent(i));
            csum += yi;

            if (sumn < 0.9 * total)
            {
                sumn += yi;
                sumx += yi * binx;
                sumxx+= yi * binx * binx;
                iend = i;
            }
        }

        const float localMean(sumx / sumn);
        const float localMeanSquared(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSquared - localMean * localMean));

        if (localRms < rmsmin)
        {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;

            if (fixDistributionCentre)
            {
                float centre = 91.2;

                if (mean > 3500.)
                {
                    throw InvalidEnergyException();
                }

                if ((mean > 2500.) && (mean < 3500.))
                    centre = 3000.;

                if ((mean > 1500.) && (mean < 2500.))
                    centre = 2000.;

                if ((mean > 700.) && (mean < 1500.))
                    centre = 1000.;

                if ((mean > 400.) && (mean < 700.))
                    centre = 500.;

                if ((mean > 250.) && (mean < 400.))
                    centre = 360.;

                if ((mean > 150.) && (mean < 250.))
                    centre = 200.;

                sigma = rms / mean * sqrt(centre);
                sigmasigma = sigma / std::sqrt(total);

                frac = rms / mean * std::sqrt(2) * 100.;
                efrac = frac / std::sqrt(total);
            }
            else
            {
                sigma = rms;
                sigmasigma = sigma / std::sqrt(total);
            }
        }
    }

    if (print)
    {
        std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries), rawrms: " << rawRms << ", rms90: " << rmsmin
                  << " (" << low << "-" << high << "), mean: " << mean << ", sigma: " << sigma << "+-" << sigmasigma;
        (fixDistributionCentre) ? (std::cout << ", sE/E: " << frac << "+-" << efrac << std::endl) : (std::cout << std::endl);
    }

    resolution = frac;
    resolutionError = efrac;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const char *InvalidEnergyException::what() const throw()
{
    return "AnalysisHelper::CalcRms - Can only fix distribution centre for small range of specific energies. ";
}

} // namespace pandora_analysis
