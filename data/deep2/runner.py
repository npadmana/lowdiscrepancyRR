import os

nrands = [1000, 4000, 16000, 64000, 256000, 1024000, 4096000, 16384000]
thetas = [0.0078125,0.015625,0.03125,0.0625,0.125,0.25]
nsim = 100
execfn = "../../cpp/build/ang2d/deep2"
nbin = len(thetas)-1


for ii in range(nbin) :
  thetamin = thetas[ii]
  thetamax = thetas[ii+1]

  for irand in nrands :
    savefn = "rrdump_%02i_%08i.dat"%(ii,irand)
    comm = "%s --nrand %i --maskfn windowf.31.fits  --thetamin %f --thetamax %f --nsim %i --save %s"%(execfn, irand, thetamin, thetamax, nsim, savefn)
    #print comm
    os.system(comm)



