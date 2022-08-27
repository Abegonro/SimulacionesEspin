/////////////////////////////////////////////////////////////////////////
//
// A simple spin analysis study (TFG Ana Belen)
// 
/////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TF1.h"

// declarations
void get_pvalue_fit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue, double& pol);
void get_pvalue_nofit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue);

// main function
void coslambdab147(Int_t numExps=50, Int_t numEvts=50000, Bool_t doFit=kTRUE, Bool_t doFitpol=kTRUE, Int_t seed=1)
{
  assert(numExps>0 && numEvts>0);
  // printout_level=0 minimum printout, !=0 print first 100 events when generating
  TRandom3 random(seed);
  printf("Generating %u experiments with %u events/experiment \n",numExps,numEvts);
  
  // Define PDFs with pol (normalized to 1)
  TF1 *pdf1 = new TF1("pdf1","[0]*0.5",-1,1); 
  TF1 *pdf2 = new TF1("pdf2","[0]*(0.5+0.25*(3.*pow(x,2)-1))",-1,1);
  TF1 *pdf3 = new TF1("pdf3","[0]*(0.5-2/7.*(3.*pow(x,2)-1)+3./112.*(3.-30.*pow(x,2)+35.*pow(x,4)))",-1,1);

  TH1F *data = new TH1F("data","data",10,-1,1);  // this makes a bin width of 0.2, ie. Entries/0.2
  double binwidth = 0.2;
  TString binwidth_str="0.2";
    
  // Loop over pseudoexperiments
  double chi2_pdf1_mean(0.),chi2_pdf2_mean(0.),chi2_pdf3_mean(0.);
  int ndof_pdf1,ndof_pdf2,ndof_pdf3;
  vector<double> pol_pdf1_list,pol_pdf2_list,pol_pdf3_list;
  double pol_pdf1_mean(0.),pol_pdf2_mean(0.),pol_pdf3_mean(0.);
  
  for (Int_t iExp=0;iExp<numExps;iExp++) {
    
    printf("Experiment %u \n",iExp);

    // clean-up
    data->Reset();

    // Generate events with pdf1 and fill in an histogram
    pdf1->SetParameter(0,1.); // set normalization to 1 
    /*
    // generate events by hand
    double x,y,ymax(pdf2->GetMaximum());
    for (Int_t iEvt=0;iEvt<numEvts;iEvt++) {
    y = ymax;
    while (y>ymax) {
    x = random.Uniform(-1.,1.);
    y = random.Uniform(0.,ymax);
    }
    data->Fill(x);
    }
    */
    // generat events automatically
    data->FillRandom("pdf1",numEvts);
    data->SetMarkerColor(kBlack);
    data->SetLineColor(kBlack);




    // Plot PDFs and generated histogram
    TCanvas* c = new TCanvas("c","c",900,900);
    //
    pdf1->SetParameter(1,0.05); // set pol to 0.05
    pdf2->SetParameter(1,0.05);
    pdf3->SetParameter(1,0.05);
    //
    pdf1->SetParameter(0,numEvts*binwidth); // set normalization to numEvts
    pdf2->SetParameter(0,numEvts*binwidth); 
    pdf3->SetParameter(0,numEvts*binwidth); 
    //
    if (!doFitpol) {
      pdf1->FixParameter(1,0.05); // fixes pol to zero
      pdf2->FixParameter(1,0.05);
      pdf3->FixParameter(1,0.05);

      pdf1->FixParameter(2, 0.05); // fixes pol to zero
      pdf2->FixParameter(2, 0.05);
      pdf3->FixParameter(2, 0.05);
    }
    {
        gStyle->SetOptStat(0);
        gStyle->SetTextFont(13);
        gStyle->SetStripDecimals(kFALSE);
        data->SetTitle("");
        data->GetXaxis()->SetTitle("cos#theta_{#Lambda_{c}}");
        data->GetYaxis()->SetTitleOffset(0);
        data->GetYaxis()->SetTitle("Número de cuentas/" + binwidth_str);
        data->GetXaxis()->SetLabelSize(0.03);
        data->GetYaxis()->SetLabelSize(0.03);
        data->SetMarkerStyle(20);
        data->SetMarkerSize(0.8);
        data->SetMarkerColor(kBlack);

        TAxis* xaxis = (TAxis*)data->GetXaxis();
        TAxis* yaxis = (TAxis*)data->GetYaxis();


        xaxis->SetLabelFont(132);
        yaxis->SetLabelFont(132);
        xaxis->SetTitleFont(132);
        yaxis->SetTitleFont(132);

        data->Draw("P0E1"); }


    { auto legend = new TLegend(0.325, 0.63, 0.675, 0.85);

    legend->AddEntry(pdf1, "Espín #frac{1}{2} #rightarrow #frac{1}{2}", "l");
    legend->SetLineColor(kBlack);
    legend->AddEntry(pdf2, "Espín #frac{1}{2} #rightarrow #frac{3}{2}"); 
    legend->AddEntry(pdf3, "Espín #frac{1}{2} #rightarrow #frac{5}{2}", "l");
    legend->AddEntry(data, "Valores de la simulación", "p");
    legend->SetTextSize(0.03);
    legend->Draw(); }

    double chi2_pdf1,pvalue_pdf1,pol_pdf1;
    pdf1->SetLineColor(kBlack);
    pdf1->SetLineWidth(1);
    if (doFit) {
      data->Fit(pdf1,"R");
      get_pvalue_fit(pdf1,data,chi2_pdf1,ndof_pdf1,pvalue_pdf1,pol_pdf1);
    } else {
      pdf1->Draw("same");
      get_pvalue_nofit(pdf1,data,chi2_pdf1,ndof_pdf1,pvalue_pdf1);
    }
    //
    double chi2_pdf2,pvalue_pdf2,pol_pdf2;
    pdf2->SetLineColor(kBlue);
    pdf2->SetLineStyle(10);
    pdf2->SetLineWidth(1);
    if (doFit) {
      data->Fit(pdf2,"R+");
      get_pvalue_fit(pdf2,data,chi2_pdf2,ndof_pdf2,pvalue_pdf2,pol_pdf2);
    } else{
      pdf2->Draw("same");
      get_pvalue_nofit(pdf2,data,chi2_pdf2,ndof_pdf2,pvalue_pdf2);
    }
    //
    double chi2_pdf3,pvalue_pdf3,pol_pdf3;
    pdf3->SetLineColor(kRed);
    pdf3->SetLineStyle(3);
    pdf3->SetLineWidth(2);
    if (doFit) {
      data->Fit(pdf3,"R+");
      get_pvalue_fit(pdf3,data,chi2_pdf3,ndof_pdf3,pvalue_pdf3,pol_pdf3);
    } else {
      pdf3->Draw("same");
      get_pvalue_nofit(pdf3,data,chi2_pdf3,ndof_pdf3,pvalue_pdf3);
    }
    //
    chi2_pdf1_mean += chi2_pdf1;
    chi2_pdf2_mean += chi2_pdf2;
    chi2_pdf3_mean += chi2_pdf3;
    //
    pol_pdf1_list.push_back(pol_pdf1);
    pol_pdf1_mean += pol_pdf1;
    pol_pdf2_list.push_back(pol_pdf2);
    pol_pdf2_mean += pol_pdf2;
    pol_pdf3_list.push_back(pol_pdf3);
    pol_pdf3_mean += pol_pdf3;    
    //
    TString file_name("analysis1_exp"); file_name += to_string(iExp); file_name += ".pdf";
    c->Print(file_name,"pdf");
    
  }
    
  // Printout summary results
  printf("\n");
  chi2_pdf1_mean /= numExps;
  chi2_pdf2_mean /= numExps;
  chi2_pdf3_mean /= numExps;
  //
  pol_pdf1_mean /= numExps;
  pol_pdf2_mean /= numExps;
  pol_pdf3_mean /= numExps;
  double pol_pdf1_variance(0.),pol_pdf2_variance(0.),pol_pdf3_variance(0.);
  for (Int_t iExp=0;iExp<numExps;iExp++) {
    pol_pdf1_variance += pow(pol_pdf1_list[iExp]-pol_pdf1_mean,2);
    pol_pdf2_variance += pow(pol_pdf2_list[iExp]-pol_pdf2_mean,2);
    pol_pdf3_variance += pow(pol_pdf3_list[iExp]-pol_pdf3_mean,2);
  }
  pol_pdf1_variance /= (numExps-1);
  pol_pdf2_variance /= (numExps-1);
  pol_pdf3_variance /= (numExps-1);
  //
  if (doFit && doFitpol)
    printf("Average pol %s: %f +- %f\n",pdf1->GetName(),pol_pdf1_mean,sqrt(pol_pdf1_variance));
  printf("Average Chi2 / ndf %s: %f / %d --> p-value: %f\n\n",pdf1->GetName(),chi2_pdf1_mean,ndof_pdf1,TMath::Prob(chi2_pdf1_mean,ndof_pdf1));
  //
  if (doFit && doFitpol)
    printf("Average pol %s: %f +- %f\n",pdf2->GetName(),pol_pdf2_mean,sqrt(pol_pdf2_variance));
  printf("Average Chi2 / ndf %s: %f / %d --> p-value: %f\n\n",pdf2->GetName(),chi2_pdf2_mean,ndof_pdf2,TMath::Prob(chi2_pdf2_mean,ndof_pdf2));
  //
  if (doFit && doFitpol)
    printf("Average pol %s: %f +- %f\n",pdf3->GetName(),pol_pdf3_mean,sqrt(pol_pdf3_variance));
  printf("Average Chi2 / ndf %s: %f / %d --> p-value: %f\n",pdf3->GetName(),chi2_pdf3_mean,ndof_pdf3,TMath::Prob(chi2_pdf3_mean,ndof_pdf3));

  return;
}


void get_pvalue_fit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue, double& pol)
{
  printf("Normalization %s: %4.3f +/- %4.3f\n",pdf->GetName(),pdf->GetParameter(0),pdf->GetParError(0));
  pol = pdf->GetParameter(1);
  printf("pol %s: %4.3f +/- %4.3f\n",pdf->GetName(),pol,pdf->GetParError(1));
  chi2 = data->GetFunction(pdf->GetName())->GetChisquare();
  ndof = data->GetFunction(pdf->GetName())->GetNDF();
  pvalue = TMath::Prob(chi2,ndof);
  printf("Chi2 / ndf %s: %f / %d --> p-value: %f\n",pdf->GetName(),chi2,ndof,pvalue);
}

void get_pvalue_nofit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue)
{
  chi2 = data->Chisquare(pdf,"L"); // Use option "R" for restricting the chisquare calculation to the given range of the function
                                   // Use option "L" for using the chisquare based on the poisson likelihood (Baker-Cousins Chisquare)
  ndof = data->GetNbinsX();
  pvalue = TMath::Prob(chi2,ndof);
  printf("Chi2 / ndf %s: %f / %d --> p-value: %f\n",pdf->GetName(),chi2,ndof,pvalue);
}
