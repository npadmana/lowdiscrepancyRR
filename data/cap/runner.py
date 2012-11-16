import os

nrands = [1000, 4000, 16000, 64000, 256000, 1024000, 4096000, 16384000]
thetamaxs = [0.1, 1.0, 10.0, 179.999]
nsim = 100
execfn = "../../cpp/build/ang2d/capintegrate"

for ii in range(len(thetamaxs)) :
  theta = thetamaxs[ii]

  for irand in nrands :
    savefn = "rrdump_%02i_%08i.dat"%(ii,irand)
    comm = "%s --nrand %i --thetamax %f --decmin 30.0 --nsim %i --save %s"%(execfn, irand, theta, nsim, savefn)
    #print comm
    os.system(comm)



