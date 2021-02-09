#!/usr/bin/python
from coffea import util,hist
import json
import os
import subprocess

import argparse

itoa = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

#python hadd_coffea.py --prefix hists_sum_bkg_ --samples 0 20 40 60 80 100 --outname hists_sum_bkg --indir Sep21_Trig -n
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', metavar='indir', default='./', help='indir')
parser.add_argument('--prefix', metavar='prefix', default='', type=str, help='prefix')
parser.add_argument('--samples', metavar='samples', help='samples', nargs='+')
parser.add_argument('--outname', metavar='outname', default='hists_sum', help='outname', type=str)
parser.add_argument('-s', '--doscale', action='store_true')
parser.add_argument('--chunk', metavar='chunk', default=5, type=int, help='chunk')
args = parser.parse_args()

indir = args.indir

chunk_size = args.chunk

onlyfiles = ["%s%s"%(args.prefix,s) for s in args.samples]

chunk_names = []
for ii,i in enumerate(range(0,len(onlyfiles),chunk_size)):
  print('Chunk',i)
  flist = []
  if (i+chunk_size<len(onlyfiles)):
    for fi in onlyfiles[i:i+chunk_size]:
      x = "%s/%s.coffea"%(indir,fi)
      flist.append(util.load(x))
  else:
    for fi in onlyfiles[i:]:
      x = "%s/%s.coffea"%(indir,fi)
      flist.append(util.load(x))

  for key in flist[0]:
    if isinstance(flist[0][key], hist.Hist):
      for fi in range(1,len(flist)):
        flist[0][key].add(flist[fi][key])
    else:
      for fi in range(1,len(flist)):
        flist[0][key] = flist[0][key] + flist[fi][key]
  
  print(flist[0])
  
  util.save(flist[0],'%s/%s_%s.coffea' % (indir,args.outname,itoa[ii]))

  for f in flist:
    del f

  chunk_names.append('%s/%s_%s.coffea' % (indir,args.outname,itoa[ii]))

print(chunk_names)

flist = [ util.load(x) for x in chunk_names ]

for key in flist[0]:
  if isinstance(key, hist.Hist):
    for fi in range(1,len(flist)):
      flist[0][key].add(flist[fi][key])
  else:
    for fi in range(1,len(flist)):
      flist[0][key] = flist[0][key] + flist[fi][key]
  
print(flist[0])

xs = {}
with open('../data/xsec.json', 'r') as f:
    xs = json.load(f)
    
scale1fb = {k: xs[k] * 1000. / w for k, w in flist[0]['sumw'].items()}
for s in args.samples:
    if s not in scale1fb: scale1fb[s] = 1.

print('doscale =',args.doscale)
if args.doscale:
    for key in flist[0]:
        if isinstance(flist[0][key], hist.Hist):
            #out[key].scale(scale(scale1fb, 'dataset'))
            flist[0][key].scale(scale1fb, 'dataset')
        else:
            print(key,flist[0][key])
            if key=='sumw':
                continue
            for samp in flist[0][key]:
                for x in flist[0][key][samp]:
                    flist[0][key][samp][x] = flist[0][key][samp][x]*scale1fb[samp]
        print(key,flist[0][key])
    
  
util.save(flist[0],'%s/%s.coffea' % (indir,args.outname))
for x in chunk_names:
  os.system("rm %s" % (indir,x))
