#!/usr/bin/python
from coffea import util,hist
import json
import os
import subprocess
from collections import defaultdict

import argparse

#python hadd_coffea.py --prefix htt_runtest_gen_zll_ --samples 200_400 400_600 600_800 800_1200 --outname htt_runtest_gen_zll
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--indir', metavar='indir', default='./', help='indir')
parser.add_argument('--prefix', metavar='prefix', default='', type=str, help='prefix')
parser.add_argument('--samples', metavar='samples', help='samples', nargs='+')
parser.add_argument('--outname', metavar='outname', default='hists_sum', help='outname', type=str)
parser.add_argument('-n', '--noscale', action='store_true')
parser.add_argument('-j', '--json', action='store_true')
args = parser.parse_args()

indir = args.indir

chunk_size = 40

onlyfiles = ["%s%s"%(args.prefix,s) for s in args.samples]

chunk_names = []
for i in range(0,len(onlyfiles),chunk_size):
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

  try: 
    for key in flist[0]:
      if isinstance(flist[0][key], hist.Hist):
        for fi in range(1,len(flist)):
          flist[0][key].add(flist[fi][key])
      else:
        for fi in range(1,len(flist)):
          flist[0][key] = flist[0][key] + flist[fi][key]
  except:
    for fi in range(1,len(flist)):
      flist[0].update(flist[fi])
  
  print(flist[0])
  
  util.save(flist[0],'%s/%s_%i.coffea' % (indir,args.outname,i))

  for f in flist:
    del f

  chunk_names.append('%s/%s_%i.coffea' % (indir,args.outname,i))

print(chunk_names)

flist = [ util.load(x) for x in chunk_names ]

try:
  for key in flist[0]:
    if isinstance(key, hist.Hist):
      for fi in range(1,len(flist)):
        flist[0][key].add(flist[fi][key])
    else:
      for fi in range(1,len(flist)):
        flist[0][key] = flist[0][key] + flist[fi][key]
except: 
  for fi in range(1,len(flist)):
    flist[0].update(flist[fi])
  
print(flist[0])

xs = {}
with open('../data/xsec.json', 'r') as f:
    xs = json.load(f)

print('noscale =',args.noscale)
if not args.noscale:
    scale1fb = {k: xs[k] * 1000. / w for k, w in flist[0]['sumw'].items()}
    for s in args.samples:
        if s not in scale1fb: scale1fb[s] = 1.
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

def makehash():
    return defaultdict(makehash)

if args.json:
    the_dict = makehash()
    for key in flist[0]:
        if isinstance(flist[0][key], hist.Hist):
            the_arr = flist[0][key].values()[()]
            the_axes = []
            for ax in flist[0][key].axes():
                the_axes.append([ax.name, ax.edges()])
            for i1 in range(len(the_axes[0][1])-1):
                for i2 in range(len(the_axes[1][1])-1):
                    the_dict[key]["ratio"]["%s:[%.2f, %.2f]"%(the_axes[0][0],the_axes[0][1][i1],the_axes[0][1][i1+1])]["%s:[%.2f, %.2f]"%(the_axes[1][0],the_axes[1][1][i2],the_axes[1][1][i2+1])]["value"] = the_arr[i1,i2]
    with open('%s/%s.json'% (indir,args.outname), 'w') as fp:
        json.dump(the_dict, fp,indent=4)
  
else:
    util.save(flist[0],'%s/%s.coffea' % (indir,args.outname))

for i,x in enumerate(chunk_names):
  os.system("rm %s/%s_%i.coffea" % (indir,args.outname,i*chunk_size))
