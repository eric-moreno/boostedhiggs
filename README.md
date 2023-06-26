# boostedhiggs

## Post-processing

Link to output is here:
```
/uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/*hist
```

### Make datacards

```
python test/makeCardsPhi.py --hist /uscms_data/d3/drankin/HTauTau/boostedhiggs_v2/condor/Nov30_2017_UL/hists_sum_ --year 2017 --lumi 41.5 --tag Dec02_2017 --label 39 --hPtBinsLep None --hPtCut 300. --hPtBinsHad None --shapeRegionsHad fail fail fail --metCutLep 75. --lowMetCutHad 75.
```
