#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;

void analysis(){

    ifstream file1 ("ARR_193nm_pol0_HPK_1.txt");
    ifstream file2 ("ARR_193nm_pol90_HPK_1.txt");
    ifstream file3 ("ARR_193nm_pol0_HPK_2.txt");
    ifstream file4 ("ARR_193nm_pol90_HPK_2.txt");
    ifstream file5 ("ARR_193nm_pol0_LF.txt");
    ifstream file6 ("ARR_193nm_pol90_LF.txt");
    ifstream file11("ARR_193nm_pol0_STD.txt");
    ifstream file12("ARR_193nm_pol90_STD.txt");
    ifstream file15("ARR_193nm_pol0_SiO2.txt");
    ifstream file16("ARR_193nm_pol90_SiO2.txt");
    string file_str_0;
    string file_str_90;
    double angle[300], reflectivity[300];
    double angle_0[300], reflectivity_0[300];
    double angle_90[300], reflectivity_90[300];

    TFile *file = new TFile ("data_193nm","recreate");
    TTree *tree = new TTree ("tree", "tree");
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
    gStyle->SetLineWidth(2);
    
    gr_file1->GetXaxis()->CenterTitle();
    gr_file1->GetYaxis()->CenterTitle();
    gr_file1->GetXaxis()->SetTitle("Angle of Incidence [#circ]");
    gr_file1->GetYaxis()->SetTitle("Specular Reflectivity [%]");
    gr_file1->GetYaxis()->SetRangeUser(0, 100);
    gr_file1->GetXaxis()->SetLimits(0, 60);
    gr_file1->SetTitle("Reflectivity VS Angle");
    gr_file1->SetLineColor(2);
    gr_file1->SetMarkerColor(2);
    gr_file2->SetLineColor(3);
    gr_file2->SetMarkerColor(3);
    gr_file3->SetLineColor(4);
    gr_file3->SetMarkerColor(4);
    gr_file6->SetLineColor(31);
    gr_file6->SetMarkerColor(31);
    gr_file8->SetLineColor(46);
    gr_file8->SetMarkerColor(46);
    TCanvas *c1 = new TCanvas("c1","Reflectivity",800,600);
    c1->cd();
    gr_file1->Draw("APL");
    gr_file2->Draw("PL");
    gr_file3->Draw("PL");
    gr_file6->Draw("PL");
    gr_file8->Draw("PL");
    TLegend *leg_0 = new TLegend (0.27, 0.83, 0.70, 0.95);
    TLegend *leg_1 = new TLegend (0.44, 0.68, 0.69, 0.84);
    TLegend *leg_2 = new TLegend (0.15, 0.68, 0.50, 0.84);
    TLegend *leg_3 = new TLegend (0.66, 0.76, 0.93, 0.84);

    leg_0->AddEntry("", "Wavelength = 193 nm","");
    leg_2->AddEntry(gr_file1, "HPK VUV4 1#","lP");
    leg_2->AddEntry(gr_file2, "HPK VUV4 2#","lP");
    leg_1->AddEntry(gr_file3, "FBK LF","lP");
    leg_1->AddEntry(gr_file6, "FBK STD","lP");
    leg_3->AddEntry(gr_file8, "FBK Wafer","lP");
    leg_0->SetTextFont(132);
    leg_1->SetTextFont(132);
    leg_2->SetTextFont(132);
    leg_3->SetTextFont(132);
    leg_0->Draw("same");
    leg_1->Draw("same");
    leg_2->Draw("same");
    leg_3->Draw("same");
    c1->Update();
    gr_file8->Write();
    c1->Write();
    tree->Fill();
    c1->Print("./pic/group_193nm.root");
    c1->Print("./pic/group_193nm.png");
    c1->Print("./pic/group_193nm.pdf");
}
