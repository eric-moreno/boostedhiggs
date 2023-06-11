# Datacards

## Make datacards

```
python datacards/make_cards.py --hist /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_ --year 2017 --cat hadhad
```

## QCD prediction

HadHad
       SR                             Inv (antilep-id): QCD_CR
    [F]    [L]    [P]              [F]    [L]    [P]
CR   A      C      D   (faildphi)   G      I      J  (qcdfaildphi)
SR   B      E      F                H      K      L  (qcdnom)

LepHad
       SR                             Inv (iso): QCD_CR
    [F]    [L]    [P]              [F]    [L]    [P]
CR   A      C      D   (lowmet)     G      I      J  (lowmet_faildphi)
SR   B      E      F                H      K      L  (qcdnom)
                                                                                                

- First method
  (Shape A_Fail_SR) * (each Pass_CR/G_Fail_CR: H/G) *  (A/A, C/A, D/A ?)
  shape * (shape region to signal region) * (correct yield in case shape is taken from [F] but need [L or P])                                                                          
  e.g.

  For HadHad SR
     hists_qcd["sig"]["faildphi"]["fail"][0] = A_Fail_SR                                                                                                                        
     qcd_ratio["sig"]["faildphi"]["fail"]["nom"] = H/G                                                                                                                                    
       temp/denom - normalized to num
       temp = hists_qcd["qcd_cr"]["nom"]["fail"] (H)
       denom = hists_qcd["qcd_cr"]["faildphi"]["fail"][0] (G)
       num = hists_qcd["sig"]["faildphi"]["fail"] (A)
     ratio_F["qcd"]["faildphi"][region] = A/A, C/A, D/A                                                                                                                                   
               
   For LepHad SR
      hists_qcd["sig"]["lowmet"]["fail"][0] = A_Fail_SR
      qcd_ratio["sig"]["lowmet"]["fail"]["nom"] = H/G                                                                                                                                    
         temp/denom - normalized to num
         temp = hists_qcd["qcd_cr"]["nom"]["fail"] (H)
         denom = hists_qcd["qcd_cr"]["lowmet"]["fail"][0] (G)
         num = hists_qcd["sig"]["lowmet"]["fail"] (A)
      ratio_F["qcd"]["lowmet"][region] = A/A, C/A, D/A

- Second method:                                                                                                                                                                       
  (Shape H_Pass_CR) * (each Fail_SR/G_Fail_CR: A/G,C/G,D/G * B/A,E/C,F/D ?) * (H/G, K/G, L/G ?)                                                                                        
  shape * (shape region to signal region) * (correct yield in case shape is taken from [F] but need [L or P]