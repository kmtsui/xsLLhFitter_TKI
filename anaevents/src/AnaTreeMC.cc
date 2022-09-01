// This is the code that actually reads int he MC tree and fills the event info.
// The tree should be produced by feeding a HL2 microtree into the treeconvert macro.

#include "AnaTreeMC.hh"

AnaTreeMC::AnaTreeMC(const std::string& file_name, const std::string& tree_name, bool extra_var)
    : read_extra_var(extra_var)
{
    fChain = new TChain(tree_name.c_str());
    fChain->Add(file_name.c_str());
    SetBranches();
}

AnaTreeMC::~AnaTreeMC()
{
    if(fChain != nullptr)
        delete fChain->GetCurrentFile();
}

long int AnaTreeMC::GetEntry(long int entry) const
{
    // Read contents of entry.
    if(fChain == nullptr)
        return -1;
    else
        return fChain->GetEntry(entry);
}

void AnaTreeMC::SetBranches()
{
    // Set branch addresses and branch pointers
    //fChain->SetBranchAddress("nutype", &nutype);
    fChain->SetBranchAddress("cutBranch", &sample);
    //fChain->SetBranchAddress("topology_multiproton", &topology);
    fChain->SetBranchAddress("topology", &topology);
    fChain->SetBranchAddress("reaction", &reaction);
    fChain->SetBranchAddress("target", &target);
    fChain->SetBranchAddress("D1True", &D1True);
    fChain->SetBranchAddress("D1Rec", &D1Reco);
    fChain->SetBranchAddress("D2True", &D2True);
    fChain->SetBranchAddress("D2Rec", &D2Reco);
    //fChain->SetBranchAddress("q2_true", &q2_true);
    //fChain->SetBranchAddress("q2_reco", &q2_reco);
    fChain->SetBranchAddress("Enutrue", &enu_true);
    //fChain->SetBranchAddress("enu_reco", &enu_reco);
    fChain->SetBranchAddress("weight", &weight);

    if(read_extra_var)
    {
        //Put extra variables here.
        fChain->SetBranchAddress("muMomRec", &muMomRec);
        fChain->SetBranchAddress("muCosThetaRec", &muCosThetaRec);
        fChain->SetBranchAddress("pMomRec", &pMomRec);
        fChain->SetBranchAddress("pCosThetaRec", &pCosThetaRec);
        fChain->SetBranchAddress("piMomRec", &piMomRec);
        fChain->SetBranchAddress("piCosThetaRec", &piCosThetaRec);

        fChain->SetBranchAddress("muMomTrue", &muMomTrue);
        fChain->SetBranchAddress("muCosThetaTrue", &muCosThetaTrue);
        fChain->SetBranchAddress("pMomTrue", &pMomTrue);
        fChain->SetBranchAddress("pCosThetaTrue", &pCosThetaTrue);
        fChain->SetBranchAddress("piMomTrue", &piMomTrue);
        fChain->SetBranchAddress("piCosThetaTrue", &piCosThetaTrue);
    }
}

void AnaTreeMC::GetEvents(std::vector<AnaSample*>& ana_samples,
                          const std::vector<SignalDef>& v_signal, const bool evt_type)
{
    if(fChain == nullptr || ana_samples.empty())
        return;

    ProgressBar pbar(60, "#");
    pbar.SetRainbow();
    pbar.SetPrefix(std::string(TAG + "Reading Events "));

    long int nentries = fChain->GetEntries();
    long int nbytes   = 0;

    std::cout << TAG << "Reading events...\n";
    for(long int jentry = 0; jentry < nentries; jentry++)
    {

        nbytes += fChain->GetEntry(jentry);

        /*if (topology<4) topology=topology; //0 CC0pi0p, 1 CC0pi1p, 2 CC0pinp, 3 CC1pi0p
        else if (topology==4&&target==1) topology=topology; //4 CC1pi1p-H
        else if (topology==4&&target==6) topology=topology+1; //5 CC1pi1p-C
        else if (topology<7) topology=topology+2; //6 CC1pi1p-other, 7 CC1pinp, 8 CCOther
        else if (topology==999) topology=9; //9 backg(NC+antinu+nue)
        else if (topology==7)   topology=10; //10 OOFV*/

        AnaEvent ev(jentry);
        ev.SetTrueEvent(evt_type);
        ev.SetFlavor(nutype);
        ev.SetSampleType(sample);
        ev.SetTopology(topology); // mectopology (i.e. CC0Pi,CC1Pi etc)
        ev.SetReaction(reaction); // reaction (i.e. CCQE,CCRES etc)
        ev.SetTarget(target);
        ev.SetTrueEnu(enu_true);
        ev.SetRecoEnu(enu_reco);
        ev.SetTrueD1(D1True);
        ev.SetRecoD1(D1Reco);
        ev.SetTrueD2(D2True);
        ev.SetRecoD2(D2Reco);
        ev.SetEvWght(weight);
        ev.SetEvWghtMC(weight);
        ev.SetQ2True(q2_true);
        ev.SetQ2Reco(q2_reco);

        if(read_extra_var)
        {
            //Put extra variables here.
            ev.SetMuMom(muMomRec);
            ev.SetMuCostheta(muCosThetaRec);
            ev.SetPrMom(pMomRec);
            ev.SetPrCostheta(pCosThetaRec);
            ev.SetPiMom(piMomRec);
            ev.SetPiCostheta(piCosThetaRec);

            ev.SetMuMomTrue(muMomTrue);
            ev.SetMuCosthetaTrue(muCosThetaTrue);
            ev.SetPrMomTrue(pMomTrue);
            ev.SetPrCosthetaTrue(pCosThetaTrue);
            ev.SetPiMomTrue(piMomTrue);
            ev.SetPiCosthetaTrue(piCosThetaTrue);
        }

        int signal_type = 0;
        for(const auto& sd : v_signal)
        {
            bool sig_passed = true;
            for(const auto& kv : sd.definition)
            {
                bool var_passed = false;
                for(const auto& val : kv.second)
                {
                    if(ev.GetEventVar(kv.first) == val)
                        var_passed = true;
                }
                sig_passed = sig_passed && var_passed;
            }
            if(sig_passed)
            {
                ev.SetSignalType(signal_type);
                ev.SetSignalEvent();
                break;
            }
            signal_type++;
        }

        for(auto& s : ana_samples)
        {
            if(s->GetSampleID() == sample)
                s->AddEvent(ev);
        }

        if(jentry % 2000 == 0 || jentry == (nentries - 1))
            pbar.Print(jentry, nentries - 1);
    }

    for(auto& sample : ana_samples)
        sample->PrintStats();
}
