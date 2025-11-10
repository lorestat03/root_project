#include <TApplication.h>
#include <TCanvas.h>
#include "Distribution.hpp"

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    Distribution dist(5.2, 1.8, 0.2, 0, 2);
    
    int nEvents = 10;
    int nBins   = 50;
    int nRuns   = 100;

    TCanvas* c1 = new TCanvas("c1", "Distribuzione media Monte Carlo", 800, 600); //distribuzione media con incertezze da rigenerazione
    dist.drawAverageHistogram(nEvents, nBins, nRuns, c1);
    TCanvas* c2 = new TCanvas("c2", "Fit distribuzione media", 800, 600);
    dist.drawAverageHistogram(10000, 50, 50, c2, false, true);
    TCanvas* c3= new TCanvas("c3", "distribuzione dopo BeanSMearing",800,600);
    dist.simulateBinSmearing(50, 1000, c3);
    TH1D* hMean = dist.getLastMeanHistogram();  
    TH1D* hSmear = dist.getLastSmearHistogram(); 
    dist.compareChi2(hMean, hSmear); //mi stampa sul terminale i 2 chi quadro per hMean e  hSmearMean
   


  app.Run();
    return 0;
}
