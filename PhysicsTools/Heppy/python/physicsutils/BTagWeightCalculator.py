import ROOT

class BTagWeightCalculator:
    def __init__(self, fn_hf, fn_lf) :
        self.pdfs = {}
        self.init(fn_hf, fn_lf)

    def init(self, fn_hf, fn_lf) :
        print "[BTagWeightCalculator]: Initializing from files", fn_hf, fn_lf

        self.pdfs["hf"] = self.getHistosFromFile(fn_hf)
        self.pdfs["lf"] = self.getHistosFromFile(fn_lf)
        return True


#    KEY: TH1D         csv_ratio_Pt2_Eta1_final
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_JESUp;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_JESDown;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_HFUp;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_HFDown;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_Stats1Up;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_Stats1Down;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_Stats2Up;1
#    KEY: TH1D         csv_ratio_Pt2_Eta1_final_Stats2Down;1

    def getHistosFromFile(self, fn):
        ret = {}
        tf = ROOT.TFile(fn)
        ROOT.gROOT.cd()
        for k in tf.GetListOfKeys():
            kn = k.GetName()
            if not kn.startswith("csv_ratio"):
                continue
            spl = kn.split("_")

            if spl[2] == "all":
                ptbin = -1
                etabin = -1
                kind = "all"
                syst = "nominal"
            else:
                ptbin = int(spl[2][2:])
                etabin = int(spl[3][3:])
                kind = spl[4]
                if len(spl)==6:
                    syst = spl[5]
                else:
                    syst = "nominal"
            ret[(ptbin, etabin, kind, syst)] = k.ReadObj().Clone()
        return ret

if __name__ == "__main__":
    bw = BTagWeightCalculator("csv/csv_rwt_hf_IT.root", "csv/csv_rwt_lf_IT.root")
#// fill the histograms (done once)
#void fillCSVhistos(TFile* fileHF, TFile* fileLF){
#
#    for( int iSys=0; iSys<9; iSys++ ){
#        for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
#        for( int iPt=0; iPt<3; iPt++ ){
#            for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
#        }
#    }
#    for( int iSys=0; iSys<5; iSys++ ){
#        for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
#    }
#
#    // CSV reweighting /// only care about the nominal ones
#    for( int iSys=0; iSys<9; iSys++ ){
#        TString syst_csv_suffix_hf = "final";
#        TString syst_csv_suffix_c = "final";
#        TString syst_csv_suffix_lf = "final";
#
#        switch( iSys ){
#        case 0:
#            // this is the nominal case
#            break;
#        case 1:
#            // JESUp
#            syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
#            syst_csv_suffix_c    = "final_cErr1Up";
#            break;
#        case 2:
#            // JESDown
#            syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
#            syst_csv_suffix_c    = "final_cErr1Down";
#            break;
#        case 3:
#            // purity up
#            syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
#            syst_csv_suffix_c    = "final_cErr2Up";
#            break;
#        case 4:
#            // purity down
#            syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
#            syst_csv_suffix_c    = "final_cErr2Down";
#            break;
#        case 5:
#            // stats1 up
#            syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
#            break;
#        case 6:
#            // stats1 down
#            syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
#            break;
#        case 7:
#            // stats2 up
#            syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
#            break;
#        case 8:
#            // stats2 down
#            syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
#            break;
#        }
#
#        for( int iPt=0; iPt<6; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );
#
#        if( iSys<5 ){
#            for( int iPt=0; iPt<6; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
#        }
#
#        for( int iPt=0; iPt<4; iPt++ ){
#            for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
#        }
#    }
#
#    return;
#}


