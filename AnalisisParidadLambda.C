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

// Declaración de funciones para obtener el p-value en el caso de ajustar o no las distribuciones de probabilidad.
void get_pvalue_nofit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue);

// Función principal (crea los pseudoexperimentos y representa el resultado en un gráfico junto con las distribuciones de probabilidad de cada espín. Devuelve el p-value para cada espín)
void AnalisisParidadLambda(Int_t numExps=10, Int_t numEvts=796, Bool_t doFit=kFALSE, Bool_t doFitBeta=kFALSE, Int_t seed=1)
{
  assert(numExps>0 && numEvts>0);
  TRandom3 random(seed);
  printf("Generating %u experiments with %u events/experiment \n",numExps,numEvts);
  
  // Definición de PDFs con beta (normalizadas a 1)
  TF1 *pdf1 = new TF1("pdf1","[0]*0.5*(1+0.45*0.61685*0.982*x)",-1,1); ////DISTRIBUCIÓN PARIDAD NEGATIVA
  TF1 *pdf2 = new TF1("pdf2","[0]*0.5*(1-0.45*0.61685*0.982*x)",-1,1); ////DISTRIBUCIÓN PARIDAD POSITIVA

  TH1F *data = new TH1F("data","data",10,-1,1);  // Define la anchura del bin para los datos generados aleatoriamente
  double binwidth = 0.2;
  TString binwidth_str="0.2";
    
  // Generación aleatoria de valores
  double chi2_pdf1_mean(0.),chi2_pdf2_mean(0.);
  int ndof_pdf1,ndof_pdf2;
 
  for (Int_t iExp=0;iExp<numExps;iExp++) {
    
    printf("Experiment %u \n",iExp);

    // clean-up
    data->Reset();
    
    // Crea los eventos con una distribución pdf2 y rellena los valores en un histograma
    pdf1->SetParameter(0,1.); // Normalización de la probabilidad a 1
    /*
    // Generar eventos a mano
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
    // Generar eventos automáticamente
    data->FillRandom("pdf1",numEvts);
    data->SetLineColor(kBlack);

    // Representar PDFs y el histograma generado
    TCanvas* c = new TCanvas("c","c",900,900);
    //
    pdf1->SetParameter(0,numEvts*binwidth); // Normaliza al número de eventos
    pdf2->SetParameter(0,numEvts*binwidth); 

   
    double chi2_pdf1,pvalue_pdf1,beta_pdf1;
    pdf1->SetLineColor(kBlack);
    pdf1->SetLineWidth(1);
    if (doFit) {
      data->Fit(pdf1,"R");
    } else {
        data->Fit(pdf1,"R");
      get_pvalue_nofit(pdf1,data,chi2_pdf1,ndof_pdf1,pvalue_pdf1);
    }
    //
    double chi2_pdf2,pvalue_pdf2,beta_pdf2;
    pdf2->SetLineColor(kBlue);
    pdf2->SetLineStyle(10);
    pdf2->SetLineWidth(1);
    if (doFit) {
      data->Fit(pdf2,"R+");
    } else{
        data->Fit(pdf2, "R+");
      get_pvalue_nofit(pdf2,data,chi2_pdf2,ndof_pdf2,pvalue_pdf2);
    }
    //
  


    // Estilo del histograma
    {
        gStyle->SetOptStat(0);
        gStyle->SetTextFont(13);
        gStyle->SetStripDecimals(kFALSE);
        data->SetTitle("");
        data->GetXaxis()->SetTitle("cos#phi_{p}");
        data->GetYaxis()->SetTitleOffset(0);
        data->GetYaxis()->SetTitle("Número de cuentas/" + binwidth_str);
        data->GetXaxis()->SetLabelSize(0.03);
        data->GetYaxis()->SetLabelSize(0.03);
        data->SetMaximum(115.);   // along          
        data->SetMinimum(35.);  //   Y    
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

    //Leyenda del histograma

    { auto legend = new TLegend(0.325, 0.63, 0.675, 0.85);

    legend->AddEntry(pdf1, "Paridad #eta=-1", "l");
    legend->SetLineColor(kBlack);
    legend->AddEntry(pdf2, "Paridad #eta=+1"); //// #frac{3}{2}
    legend->AddEntry(data, "Valores de la simulación", "p");
    legend->SetTextSize(0.03);
    legend->Draw(); }

    //
    chi2_pdf1_mean += chi2_pdf1;
    chi2_pdf2_mean += chi2_pdf2;
    //
    TString file_name("analysisOmegabetalibre_exp"); file_name += to_string(iExp); file_name += ".pdf";
    c->Print(file_name,"pdf");
    
  }
    
  // Muestra en pantalla un resumen de los resultados
  printf("\n");
  chi2_pdf1_mean /= numExps;
  chi2_pdf2_mean /= numExps;


  printf("Average Chi2 / ndf %s: %f / %d --> p-value: %f\n\n",pdf1->GetName(),chi2_pdf1_mean,ndof_pdf1,TMath::Prob(chi2_pdf1_mean,ndof_pdf1));
  //

  printf("Average Chi2 / ndf %s: %f / %d --> p-value: %f\n\n",pdf2->GetName(),chi2_pdf2_mean,ndof_pdf2,TMath::Prob(chi2_pdf2_mean,ndof_pdf2));
  //
  return;
}


void get_pvalue_nofit(TF1 *pdf, TH1F *data, double& chi2, int& ndof, double& pvalue)
{
  chi2 = data->Chisquare(pdf,"L"); // Use option "R" for restricting the chisquare calculation to the given range of the function
                                   // Use option "L" for using the chisquare based on the poisson likelihood (Baker-Cousins Chisquare)
  ndof = data->GetNbinsX();
  pvalue = TMath::Prob(chi2,ndof);
  printf("Chi2 / ndf %s: %f / %d --> p-value: %f\n",pdf->GetName(),chi2,ndof,pvalue);
}
