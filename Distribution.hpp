#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP
#include <TF1.h>
#include <vector>
#include <TCanvas.h>
#include <TH1D.h>



class Distribution {
private:
    double k;
    double phi;
    double b;
    double xmin;
    double xmax;
    TF1* func;
    mutable TH1D* lastMean = nullptr;
    mutable TH1D* lastSmear = nullptr;

public:
    Distribution(double k_=5.2, double phi_=1.8, double b_=0.2, double xmin_=0, double xmax_=2);


       // Genera nRun istogrammi Monte Carlo e ritorna i contenuti dei bin
    std::vector<std::vector<double>> generateMultipleRuns(int nEvents, int nBins, int nRuns) const;

    // Ritorna i valori teorici della funzione nei centri dei bin
    std::vector<double> getFunctionValuesInBins(int nBins) const;
    // Calcola la media dei bin su più run
std::vector<double> computeBinMeans(const std::vector<std::vector<double>>& allBinContents) const;

// Calcola la deviazione standard dei bin su più run
std::vector<double> computeBinStdDev(const std::vector<std::vector<double>>& allBinContents, const std::vector<double>& binMeans) const;



void drawAverageHistogram(int nEvents, int nBins, int nRuns, TCanvas* canvas,bool drawMean=true, bool doFit=false) const;

void simulateBinSmearing(int nBins, int nSmear,TCanvas* canvas) const;
double compareChi2(TH1D* h1, TH1D* h2) const;
TH1D* getLastMeanHistogram() const { return lastMean; } // lastMean da salvare dentro drawAverageHistogram
TH1D* getLastSmearHistogram() const { return lastSmear; } // idem

};



#endif
