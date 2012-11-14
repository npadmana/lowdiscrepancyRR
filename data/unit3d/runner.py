import os

nrands = [1000, 4000, 16000, 64000, 256000, 1024000, 4096000, 16384000]
rbins = [0.015625, 0.03125, 0.0625, 0.125, 0.25, 0.5]
nr = len(rbins)-1
nsim = 100
execfn = "../../cpp/build/unit3d/unit3d"

for ii in range(nr) :
  rmin = rbins[ii]
  rmax = rbins[ii+1]

  for irand in nrands :
    savefn = "rrdump_%02i_%08i.dat"%(ii,irand)
    comm = "%s --nrand %i --rmin %f --rmax %f --nsim %i --save %s"%(execfn, irand, rmin, rmax, nsim, savefn)
    print comm
    os.system(comm)
    print "-----------------------------"



