# HN3 and HN5 Phase diagram
# This code is to an example to use HNRG_module.py
# It may take 5s ~ 16s to run

import numpy as np
import matplotlib.pyplot as plt
import HNRG as hnrg

def main():
  nt = 100
  ms = np.linspace(0.01, 0.64, num = nt)
  n = 5000
  l = 64
  ys = np.array([0, 0.5, 1, 2])
  y = 0
  
  fsize = 12 
  msize = 1
  
  plt.figure(1, figsize = (10,10))
  hp = hnrg.hp6() # use of class hp6()
  for yi in range(ys.size):
    y = ys[yi]
    plt.subplot(int(np.sqrt(ys.size)), int(ys.size/np.sqrt(ys.size)), yi+1)
    for i in range(ms.size):
      m =1/ ms[i]
      ks = np.unique(hp.rg_n_k(m, y, n, l))
      #plt.plot(ms[i]*np.ones((ks.size,)), ks, 'ks', markersize = 1)
      plt.semilogy(ms[i]*np.ones((ks.size,)), ks, 'ks', markersize = msize)
    plt.grid('on')
    #plt.ylim([0.01, 100])
    plt.xlim([np.amin(ms), np.amax(ms)])
    #plt.legend(loc = 2)
    plt.ylabel('$\kappa$', fontsize = fsize)
    if yi==ys.size-1:
      plt.xlabel('$1/\mu$', fontsize = fsize)
    plt.title('HP6,y='+str(y)+',steps:'+str(n-l)+'~'+str(n), fontsize = fsize)
  #plt.savefig('HNRG_eg_'+str(n)+'l'+str(l)+'.pdf')
  plt.show()
  
if __name__ == '__main__':
  main()
  

