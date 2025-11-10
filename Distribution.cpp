#include "Distribution.hpp"
#include <TF1.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TMath.h>
#include <cmath>
#include <vector>
#include <TLatex.h>
#include <TLegend.h>


Distribution::Distribution(double k_, double phi_, double b_, double xmin_, double xmax_)
    : k(k_), phi(phi_), b(b_), xmin(xmin_), xmax(xmax_)
{
    // Funzione analitica ROOT
    func = new TF1("func", "cos(5.2*x + 1.8)^2 + 0.2", xmin, xmax);
}


// Genera N run di eventi distribuiti secondo la funzione TF1

std::vector<std::vector<double>> Distribution::generateMultipleRuns(int nEvents, int nBins, int nRuns) const {
    std::vector<std::vector<double>> allBinContents(nRuns, std::vector<double>(nBins, 0.0));
    double fmin = xmin;
    double fmax = xmax;

    for (int run = 0; run < nRuns; ++run) {
        TH1D hist("hist_tmp", "Distribuzione temporanea", nBins, xmin, xmax);

        for (int i = 0; i < nEvents; ++i) {
            double x = func->GetRandom(xmin, xmax);  // ROOT genera secondo f(x)
            hist.Fill(x);
        }

        // Normalizzazione all’area unitaria (densità di probabilità)
        hist.Scale(1.0 / hist.Integral("width"));

        for (int b = 1; b <= nBins; ++b)
            allBinContents[run][b - 1] = hist.GetBinContent(b);
    }

    return allBinContents;
}


// Valori teorici della funzione analitica ai centri dei bin

std::vector<double> Distribution::getFunctionValuesInBins(int nBins) const {
    std::vector<double> values(nBins, 0.0);
    double binWidth = (xmax - xmin) / nBins;

    for (int i = 0; i < nBins; ++i) {
        double xCenter = xmin + (i + 0.5) * binWidth;
        values[i] = func->Eval(xCenter);
    }

    // Normalizziamo la funzione analitica all’integrale unitario
    double norm = func->Integral(xmin, xmax);
    for (auto &v : values) v /= norm;

    return values;
}


// Calcola le medie dei bin su più run

std::vector<double> Distribution::computeBinMeans(const std::vector<std::vector<double>>& allBinContents) const {
    int nRuns = allBinContents.size();
    int nBins = allBinContents[0].size();
    std::vector<double> means(nBins, 0.0);

    for (int b = 0; b < nBins; ++b) {
        std::vector<double> binValues(nRuns);
        for (int r = 0; r < nRuns; ++r)
            binValues[r] = allBinContents[r][b];
        means[b] = TMath::Mean(nRuns, binValues.data());
    }
    return means;
}


// Calcola la deviazione standard dei bin su più run

std::vector<double> Distribution::computeBinStdDev(const std::vector<std::vector<double>>& allBinContents,
                                                  const std::vector<double>& binMeans) const {
    int nRuns = allBinContents.size();
    int nBins = allBinContents[0].size();
    std::vector<double> stddev(nBins, 0.0);

    for (int b = 0; b < nBins; ++b) {
        std::vector<double> binValues(nRuns);
        for (int r = 0; r < nRuns; ++r)
            binValues[r] = allBinContents[r][b];
        stddev[b] = TMath::StdDev(nRuns, binValues.data());
    }
    return stddev;}
    void Distribution::drawAverageHistogram(int nEvents, int nBins, int nRuns, TCanvas* canvas,bool drawMean, bool doFit) const {
    //  Disegno la funzione analitica
    func->SetLineColor(kBlue);
    func->SetLineWidth(2);
    func->SetTitle("Distribuzione media Monte Carlo; x; f(x)");
    func->Draw();

    // Genero le run Monte Carlo
    auto allBinContents = generateMultipleRuns(nEvents, nBins, nRuns);
    auto binMeans  = computeBinMeans(allBinContents);
    auto binStdDev = computeBinStdDev(allBinContents, binMeans);

    // Creo l'istogramma medio
    TH1D* hMean = new TH1D("hMean", "Distribuzione Monte Carlo Media", nBins, xmin, xmax);
    lastMean = hMean;

    for (int i = 0; i < nBins; ++i) {
        hMean->SetBinContent(i + 1, binMeans[i]);
        hMean->SetBinError(i + 1, binStdDev[i]);
    }

    // Normalizzazione (opzionale)
    double integralFunc = func->Integral(xmin, xmax);
    double integralHist = hMean->Integral("width");
    if (integralHist > 0)
        hMean->Scale(integralFunc / integralHist);

    // Stile e disegno
    hMean->SetLineColor(kRed);
    hMean->SetLineWidth(2);
    hMean->SetMarkerStyle(20);
    hMean->SetMarkerColor(kRed);
    hMean->Draw("HIST SAME");
    hMean->Draw("E SAME");

double meanStdDev = 0;
    for (double s : binStdDev) meanStdDev += s;
    meanStdDev /= binStdDev.size();
     double meanStdDevRMS = 0;
for (double s : binStdDev) meanStdDevRMS += s*s;
meanStdDevRMS = std::sqrt(meanStdDevRMS / binStdDev.size());
    

    // Scrittura dell'incertezza media sulla canvas 
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.15, 0.85, Form("Incertezza media da rigenerazione = %.3g", meanStdDev));
    latex.DrawLatex(0.15, 0.8, Form("Incertezza RMS = %.3g", meanStdDevRMS));

    // Legenda 
    
    canvas->SaveAs("distribuzione_media.png");
     
    
    if (doFit) {
        TF1* fitFunc = new TF1("fitFunc", "[0]*cos([1]*x + [2])^2 + [3]", xmin, xmax);
        fitFunc->SetParameters(1.0, 5.0, 1.5, 0.2);
        hMean->Fit(fitFunc, "R");

   

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf = %.2f / %d",
        fitFunc->GetChisquare(), fitFunc->GetNDF()));

        for (int i = 0; i < fitFunc->GetNpar(); ++i)
        latex.DrawLatex(0.15, 0.8 - 0.05 * i,
        Form("p_{%d} = %.3f #pm %.3f",
        i, fitFunc->GetParameter(i), fitFunc->GetParError(i)));

        TLegend* legendFit = new TLegend(0.55, 0.7, 0.88, 0.88);
        legendFit->AddEntry(fitFunc, "Fit libero", "l");
        legendFit->Draw();
       hMean->Draw("E1");
     fitFunc->SetLineColor(kGreen + 2);
        fitFunc->SetLineWidth(2);
        fitFunc->Draw( "SAME");
    }

    canvas->Update();
}

void Distribution::simulateBinSmearing(int nBins, int nSmear,TCanvas* canvas ) const{
    // Creo l'istogramma della funzione teorica
    TH1D* hTrue = new TH1D("hTrue", "Bin Smearing - Distribuzione teorica e fluttuata; x; f(x)", nBins, xmin, xmax);
    lastSmear = hTrue;

    double binWidth = (xmax - xmin) / nBins;

    for (int i = 1; i <= nBins; ++i) {
        double xCenter = xmin + (i - 0.5) * binWidth;
        double val = func->Eval(xCenter);
        hTrue->SetBinContent(i, val);
    }

    // Normalizzo la funzione
    hTrue->Scale(1.0 / hTrue->Integral("width"));

    // Creo un istogramma per le fluttuazioni (media su più smearing)
    TH1D* hSmearMean = (TH1D*)hTrue->Clone("hSmearMean");
    hSmearMean->Reset();

    TRandom3 rnd(0); // generatore casuale

    // Vettore per accumulare valori fluttuati per ogni bin
    std::vector<std::vector<double>> smearedValues(nBins, std::vector<double>(nSmear, 0.0));

    for (int s = 0; s < nSmear; ++s) {
        for (int b = 1; b <= nBins; ++b) {
            double mean = hTrue->GetBinContent(b);
            double sigma = 0.05 * mean; // fluttuazione del 5% (puoi modificare)
            double fluctuated = rnd.Gaus(mean, sigma);
            smearedValues[b - 1][s] = fluctuated;
        }
    }

    // Calcolo media e deviazione standard per ogni bin
    std::vector<double> meanVals(nBins), stdVals(nBins);
    for (int b = 0; b < nBins; ++b) {
        meanVals[b] = TMath::Mean(nSmear, smearedValues[b].data());
        stdVals[b]  = TMath::StdDev(nSmear, smearedValues[b].data());
        hSmearMean->SetBinContent(b + 1, meanVals[b]);
        hSmearMean->SetBinError(b + 1, stdVals[b]);
    }

    // Disegno
    canvas->cd();
    hTrue->SetLineColor(kBlue);
    hTrue->SetLineWidth(2);
    hTrue->Draw("HIST");

    hSmearMean->SetMarkerStyle(20);
    hSmearMean->SetMarkerColor(kRed);
    hSmearMean->Draw("E1 SAME");

    // Legenda
    TLegend* leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hTrue, "Distribuzione teorica", "l");
    leg->AddEntry(hSmearMean, "Distribuzione fluttuata (Bin Smearing)", "lep");
    leg->Draw();

    // Stima incertezza media
    double meanSigma = 0;
    for (double s : stdVals) meanSigma += s;
    meanSigma /= stdVals.size();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.85, Form("Incertezza media (smearing) = %.3g", meanSigma));

    canvas->SaveAs("bin_smearing.png");
}

 double Distribution::compareChi2(TH1D* h1, TH1D* h2) const {
     double chi2ROOT = h1->Chi2Test(h2, "CHI2");
    int nBins = h1->GetNbinsX();
    double chi2Manual = 0;
    for (int i = 1; i <= nBins; ++i) {
        double diff = h1->GetBinContent(i) - h2->GetBinContent(i);
        double err = h1->GetBinError(i);
        if (err > 0) chi2Manual += diff*diff / (err*err);
    }
    std::cout << "Chi² ROOT = " << chi2ROOT << std::endl;
    std::cout << "Chi² manuale = " << chi2Manual << std::endl;
    return chi2Manual;
}
