#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

void analysis(){

    ifstream file1("Reflectance_HPK1#_pol0_57.txt");
    ifstream file2("Reflectance_HPK1#_pol90_57.txt");
    ifstream file3("Reflectance_HPK2#_pol0_57.txt");
    ifstream file4("Reflectance_HPK2#_pol90_57.txt");
    ifstream file5("Reflectance_LF_pol0_57.txt");
    ifstream file6("Reflectance_LF_pol90_57.txt");
    ifstream file11("Reflectance_STD_pol0_57.txt");
    ifstream file12("Reflectance_STD_pol90_57.txt");
    ifstream file15("Reflectance_Wafer_pol0_57.txt");
    ifstream file16("Reflectance_Wafer_pol90_57.txt");
   string file_str_0;
    string file_str_90;
    double angle[500], reflectivity[500];
    double angle_0[500], reflectivity_0[500];
    double angle_90[500], reflectivity_90[500];
    
    int npointfile=0;
    while (getline(file1, file_str_0) && getline(file2, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file1 = new TGraph(npointfile, angle, reflectivity);
  
    npointfile=0;
    while (getline(file3, file_str_0) && getline(file4, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file2 = new TGraph(npointfile, angle, reflectivity);
      
    npointfile=0;
    while (getline(file5, file_str_0) && getline(file6, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file3 = new TGraph(npointfile, angle, reflectivity);
  /* 
    npointfile=0;
    while (getline(file7, file_str_0) && getline(file8, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file4 = new TGraph(npointfile, angle, reflectivity);
   
    npointfile=0;
    while (getline(file9, file_str_0) && getline(file10, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file5 = new TGraph(npointfile, angle, reflectivity);
 */
    npointfile=0;
    while (getline(file11, file_str_0) && getline(file12, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file6 = new TGraph(npointfile, angle, reflectivity);
/*
    npointfile=0;
    while (getline(file13, file_str_0) && getline(file14, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file7 = new TGraph(npointfile, angle, reflectivity);
*/
    npointfile=0;
    while (getline(file15, file_str_0) && getline(file16, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file8 = new TGraph(npointfile, angle, reflectivity);
/*
    npointfile=0;
    while (getline(file17, file_str_0) && getline(file18, file_str_90)){
        istringstream ss_0(file_str_0);
        istringstream ss_90(file_str_90);
        ss_0>>angle_0[npointfile]>>reflectivity_0[npointfile];
        ss_90>>angle_90[npointfile]>>reflectivity_90[npointfile];
        angle[npointfile] = (angle_0[npointfile] + angle_90[npointfile])/2;
        reflectivity[npointfile] = (reflectivity_0[npointfile] + reflectivity_90[npointfile])/2;
        npointfile++;
    }
    TGraph *gr_file9 = new TGraph(npointfile, angle, reflectivity);
  
    npointfile=0;
    while (getline(file11, file_str)){
        istringstream ss(file_str);
        ss>>angle[npointfile]>>reflectivity[npointfile];
        npointfile++;
    }
    TGraph *gr_file11 = new TGraph(npointfile, angle, reflectivity);
  */  
    gStyle->SetLineWidth(2);
    
    gr_file1->GetXaxis()->CenterTitle();
    gr_file1->GetYaxis()->CenterTitle();
    gr_file1->GetXaxis()->SetTitle("Wavelength [nm]");
    gr_file1->GetYaxis()->SetTitle("Specular Reflectivity [%]");
    gr_file1->GetYaxis()->SetRangeUser(0, 90);
    gr_file1->GetXaxis()->SetLimits(110, 310);
    gr_file1->SetTitle("Reflectivity VS Wavelength");
    gr_file1->SetLineColor(2);
    //gr_file1->SetMarkerColor(2);
    gr_file2->SetLineColor(3);
    //gr_file2->SetMarkerColor(3);
    gr_file3->SetLineColor(4);
    //gr_file3->SetMarkerColor(4);
/*    gr_file4->SetLineColor(1);
    gr_file4->SetMarkerColor(1);
    gr_file5->SetLineColor(7);
    gr_file5->SetMarkerColor(7);
*/    gr_file6->SetLineColor(31);
    //gr_file6->SetMarkerColor(31);
  //  gr_file7->SetLineColor(43);
  //  gr_file7->SetMarkerColor(43);
    gr_file8->SetLineColor(46);
    //gr_file8->SetMarkerColor(46);
/*  gr_file10->SetLineColor(38);
    gr_file10->SetMarkerColor(38);
    gr_file11->SetLineColor(30);
    gr_file11->SetMarkerColor(30);
*/  
    TCanvas *c1 = new TCanvas("c1","Reflectivity",800,600);
    c1->cd();
    gr_file1->Draw("AL");
    gr_file2->Draw("L");
    gr_file3->Draw("L");
    //gr_file4->Draw("PL");
    //gr_file5->Draw("PL");
    gr_file6->Draw("L");
    //gr_file7->Draw("PL");
    gr_file8->Draw("L");
/*    gr_file10->Draw("PL");
    gr_file11->Draw("PL");
*/
    TLegend *leg_0 = new TLegend (0.27, 0.83, 0.70, 0.95);
    TLegend *leg_1 = new TLegend (0.44, 0.68, 0.69, 0.84);
    TLegend *leg_2 = new TLegend (0.15, 0.68, 0.50, 0.84);
    TLegend *leg_3 = new TLegend (0.66, 0.76, 0.93, 0.84);

    leg_0->AddEntry("", "Incident Angle = 57 #circ","");
    //leg->AddEntry("", "Pixel Size = 75 um","");
    leg_2->AddEntry(gr_file1, "HPK VUV4 1#","l");
    leg_2->AddEntry(gr_file2, "HPK VUV4 2#","l");
    leg_1->AddEntry(gr_file3, "FBK LF","l");
    //leg->AddEntry(gr_file4, "FBK RGB1#","l");
    //leg->AddEntry(gr_file5, "FBK RGB2#","l");
    leg_1->AddEntry(gr_file6, "FBK STD","l");
    //leg->AddEntry(gr_file7, "FBK Wafer Si_{3}N_{4}","l");
    leg_3->AddEntry(gr_file8, "FBK Wafer","l");
 /* leg->AddEntry(gr_file10, "Wavelength = 200 nm","l");
    leg->AddEntry(gr_file11, "Wavelength = 250 nm","l");
  */  leg_0->SetTextFont(132);
    leg_1->SetTextFont(132);
    leg_2->SetTextFont(132);
    leg_3->SetTextFont(132);
    leg_0->Draw("same");
    leg_1->Draw("same");
    leg_2->Draw("same");
    leg_3->Draw("same");
    c1->Update();
    c1->Print("./pic/group_57deg_new.root");
    c1->Print("./pic/group_57deg_new.png");
    c1->Print("./pic/group_57deg_new.pdf");
}
