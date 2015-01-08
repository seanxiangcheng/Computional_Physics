"""
Renormalization Group Recurisive Equations of Hierarchical Networks;
This module is for Xiang Cheng's dissertation research at Emory University;
For more info about this topic: 
Renormalization group: http://journals.aps.org/pre/abstract/10.1103/PhysRevE.83.021103
Hierarchical Networks (HN3, HN5, HNNP, HN6):http://iopscience.iop.org/1751-8121/41/33/335003 

Applicable Complex Networks of this module: HN3, HN5, HNNP, HN6;
Classes: 
  h3: HN3; (h3 call functions in h35)
  h5: HN5; (h5 call functions in h35)
  h35: HN3 and HN5; (just for the convenience of coding)
  hp: HNNP; (hp call functions in hp6)
  h6: HN6;  (h6 call functions in hp6)
  hp6: HNNP and HN6; (just for the convenience of coding)
  
Functions: See the definitions of each functions

Comments: This module basically consists of all the classes, functions, 
          varialbes for numerical calculations;
          Plotting and analysis are conducted in other scripts;
          Example of the usage of this module is "hnrg_eg.py";
          You are welcome to modify the code for your purpose.
          The author assumes no responsibility.
"""
import numpy as np

#*******************   Classes   ******************#
# class hp6 for Networks: HNNP and HN6
class hp6(object):
  # initialization function
  def __init__(self):
    return ;
  # recursive equation of kappa
  def kap(self, m, k, lam, y):
    value = k*lam*((1.+m)**2)/((1.+m*k)**2)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k, y):
    value = m**(2.*y)*((m+k)**2)/((1.+m*k)**2)
    return value
  
  # recursive equation of C (constant parameter)
  def con(self, m, k, c, y):
    value = c**2 * k*m**2/(1.+m)**2/(m+k)/(1.+m*k)
    return value
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, y, n):
    # initial values
    m = float(m); n = int(n) # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam 
    if n==0:
      return k
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return k
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, y, n, l):
    # initial values
    m = float(m); n = int(n); l = int(l); # to avoid data type errors
    if l>n:
      l = n
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    ks = np.zeros((l, ))
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
      if i > (n-l-1):
        ks[i-n+l] = k
    return ks
  
  # find all the k through our the iterations
  def rg_all_k(self, m, y, n):
    # initial values
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    ks = np.zeros((n+1, ))
    # RG iterations
    for i in range(n+1):
      ks[i] = k
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return ks
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, y, n, l):
    m = float(m); n = int(n); # to avoid data type errors
    if l>n:
      l = n
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    kls = np.zeros((l, 2))
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
      if i > (n-l-1):
        kls[i-n+l, 0] = k
        kls[i-n+l, 1] = lam
    return kls

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, y, n):
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    kls = np.zeros((n+1, 2))
    # RG iterations
    for i in range(n+1):
      kls[i, 0] = k
      kls[i, 1] = lam
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return kls
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, y, n):
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    c = 0
    #con(self, m, k, c, y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    c0 = c
    values = np.zeros((n+1, 3))
    # RG iterations
    for i in range(n+1):
      values[i, 0] = k
      values[i, 1] = lam
      values[i, 3] = c
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      c = self.con(m, k0, c0, y)
      k0 = k
      lam0 = lam
      c0 = c
    return values
  
  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, y, n, l):
    mm = np.zeros((2,))
    ks = self.rg_n_k(m, y, n, l)
    mm[0] = np.amin(ks)
    mm[1] = np.amax(ks)
    return mm

# class h35 for Networks: HN3 and HN5    
class h35(object):
  # initialization function
  def __init__(self):
    return ;
  # recursive equation of kappa
  def kap(self, m, k, lam, y):
    value = k*lam*2.0*(1.0+m)/(1.0+2.0*m*k+k**2)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k, y):
    value = m**(2.0*y)*((1.+k)**2*(1.+m))/(2*(1.+2.*m*k+k**2))
    return value
  
  # recursive equation of C (constant parameter)
  def con(self, m, k, c, y):
    value = c**2 * k*m/np.sqrt(2*(1.+2.*m*k+k**2)*(1.+m)**3)/(1.+k)
    return value
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, y, n):
    # initial values
    m = float(m); n = int(n) # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam 
    if n==1:
      return k
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return k
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, y, n, l):
    # initial values
    m = float(m); n = int(n); l = int(l); # to avoid data type errors
    if l>n:
      l = n
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    ks = np.zeros((l, ))
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
      if i > (n-l-1):
        ks[i-n+l] = k
    return ks
  
  # find all the k through our the iterations
  def rg_all_k(self, m, y, n):
    # initial values
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    ks = np.zeros((n+1, ))
    # RG iterations
    for i in range(n+1):
      ks[i] = k
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return ks
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, y, n, l):
    m = float(m); n = int(n); # to avoid data type errors
    if l>n:
      l = n
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    kls = np.zeros((l, 2))
    # RG iterations
    for i in range(n):
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
      if i > (n-l-1):
        kls[i-n+l, 0] = k
        kls[i-n+l, 1] = lam
    return kls

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, y, n):
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    kls = np.zeros((n+1, 2))
    # RG iterations
    for i in range(n+1):
      kls[i, 0] = k
      kls[i, 1] = lam
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      k0 = k
      lam0 = lam
    return kls
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, y, n):
    m = float(m); n = int(n); # to avoid data type errors
    k = m**2
    lam = m**(2*y)
    c = 0
    #con(self, m, k, c, y)
    # values to store old values after RG
    k0 = k
    lam0 = lam
    c0 = c
    values = np.zeros((n+1, 3))
    # RG iterations
    for i in range(n+1):
      values[i, 0] = k
      values[i, 1] = lam
      values[i, 3] = c
      lam = self.lam(m, k0, y)
      k = self.kap(m, k0, lam0, y)
      c = self.con(m, k0, c0, y)
      k0 = k
      lam0 = lam
      c0 = c
    return values
  
  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, y, n, l):
    mm = np.zeros((2,))
    ks = self.rg_n_k(m, y, n, l)
    mm[0] = np.amin(ks)
    mm[1] = np.amax(ks)
    return mm
    
# class hp for network: HNNP; 
# This class is just for convenience; 
# It recalls functions in class hp6 with y = 0
class hp(object):
    # initialization function
  def __init__(self):
    return ;
    
  # recursive equation of kappa
  def kap(self, m, k, lam):
    y = 0
    tem = hp6()
    value = tem.kap(m, k, lam, y)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k):
    y = 0
    tem = hp6()
    value = tem.lam(m, k, y)
    return value 
    
  # recursive equation of C (constant parameter)
  def con(self, m, k, c):
    y = 0
    tem = hp6()
    value = tem.con(m, k, c, y)
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, n):
    y = 0
    tem = hp6()
    return tem.rg_k(m, y, n)
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, n, l):
    y = 0
    tem = hp6()
    return tem.rg_n_k(m, y, n, l)
  
  # find all the k through our the iterations
  def rg_all_k(self, m, n):
    y = 0
    tem = hp6()
    return tem.rg_all_k(m, y, n)
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, n, l):
    y = 0
    tem = hp6()
    return tem.rg_n_kl(m, y, n, l)

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, n):
    y = 0
    tem = hp6()
    return tem.rg_all_kl(m, y, n)
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, n):
    y = 0
    tem = hp6()
    return tem.rg_all(m, y, n)
  
  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, n, l):
    y = 0
    tem = hp6()
    return tem.rg_mm_k(m, y, n, l)


# class hp for Network: HN6; 
# This class is just for convenience; 
# It recalls functions in class hp6 with y = 0
class h6(object):
    # initialization function
  def __init__(self):
    return ;
    
  # recursive equation of kappa
  def kap(self, m, k, lam):
    y = 1
    tem = hp6()
    value = tem.kap(m, k, lam, y)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k):
    y = 1
    tem = hp6()
    value = tem.lam(m, k, y)
    return value 
    
  # recursive equation of C (constant parameter)
  def con(self, m, k, c):
    y = 1
    tem = hp6()
    value = tem.con(m, k, c, y)
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, n):
    y = 1
    tem = hp6()
    return tem.rg_k(m, y, n)
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, n, l):
    y = 1
    tem = hp6()
    return tem.rg_n_k(m, y, n, l)
  
  # find all the k through our the iterations
  def rg_all_k(self, m, n):
    y = 1
    tem = hp6()
    return tem.rg_all_k(m, y, n)
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, n, l):
    y = 1
    tem = hp6()
    return tem.rg_n_kl(m, y, n, l)

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, n):
    y = 1
    tem = hp6()
    return tem.rg_all_kl(m, y, n)
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, n):
    y = 1
    tem = hp6()
    return tem.rg_all(m, y, n)
  
  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, n, l):
    y = 1
    tem = hp6()
    return tem.rg_mm_k(m, y, n, l)

# class h3 for Network: HN3; 
# This class is just for convenience; 
# It recalls functions in class h35 with y = 0
class h3(object):
    # initialization function
  def __init__(self):
    return ;
    
  # recursive equation of kappa
  def kap(self, m, k, lam):
    y = 0
    tem = h35()
    value = tem.kap(m, k, lam, y)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k):
    y = 0
    tem = h35()
    value = tem.lam(m, k, y)
    return value 
    
  # recursive equation of C (constant parameter)
  def con(self, m, k, c):
    y = 0
    tem = h35()
    value = tem.con(m, k, c, y)
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, n):
    y = 0
    tem = h35()
    return tem.rg_k(m, y, n)
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, n, l):
    y = 0
    tem = h35()
    return tem.rg_n_k(m, y, n, l)
  
  # find all the k through our the iterations
  def rg_all_k(self, m, n):
    y = 0
    tem = h35()
    return tem.rg_all_k(m, y, n)
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, n, l):
    y = 0
    tem = h35()
    return tem.rg_n_kl(m, y, n, l)

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, n):
    y = 0
    tem = h35()
    return tem.rg_all_kl(m, y, n)
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, n):
    y = 0
    tem = h35()
    return tem.rg_all(m, y, n)

  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, n, l):
    y = 0
    tem = h35()
    return tem.rg_mm_k(m, y, n, l)


# class hp for HN5; For convenience, recall functions in h35 with y = 1
class h5(object):
    # initialization function
  def __init__(self):
    return ;
    
  # recursive equation of kappa
  def kap(self, m, k, lam):
    y = 1
    tem = h35()
    value = tem.kap(m, k, lam, y)
    return value 
  
  # recursive equation of lambda
  def lam(self, m, k):
    y = 1
    tem = h35()
    value = tem.lam(m, k, y)
    return value 
    
  # recursive equation of C (constant parameter)
  def con(self, m, k, c):
    y = 1
    tem = h35()
    value = tem.con(m, k, c, y)
  
  # find the last value of kappa after n RG steps 
  def rg_k(self, m, n):
    y = 1
    tem = h35()
    return tem.rg_k(m, y, n)
    
  # find the last value of kappa after n RG steps 
  def rg_n_k(self, m, n, l):
    y = 1
    tem = h35()
    return tem.rg_n_k(m, y, n, l)
  
  # find all the k through our the iterations
  def rg_all_k(self, m, n):
    y = 1
    tem = h35()
    return tem.rg_all_k(m, y, n)
    
  # find the last l steps: kappa and lambda
  def rg_n_kl(self, m, n, l):
    y = 1
    tem = h35()
    return tem.rg_n_kl(m, y, n, l)

  # find all the kappa and lambda including the intial values
  def rg_all_kl(self, m, n):
    y = 1
    tem = h35()
    return tem.rg_all_kl(m, y, n)
  
  # find all the kappa, lambda, C, which is usually not used.
  def rg_all(self, m, n):
    y = 1
    tem = h35()
    return tem.rg_all(m, y, n)
  # find the max and min of kappa of the last l steps among n
  def rg_mm_k(self, m, n, l):
    y = 1
    tem = h35()
    return tem.rg_mm_k(m, y, n, l)


#********************** Functions to serve the classes ***********************#

# Find the histogram of the array of Fixed Points
# density True returns the normalized histogram (0, 1); 
#         False returns number of visits; 
# nbins: number of bins
# mm is the array of pre-defined [max, min] of fixed points;
# it can help to returns histogram with the same bins.
def hist_fx(fps, nbins = 20, density = True, mm = np.array([])):
  if fps.ndim > 1:
    print 'Array of fixed point has to be 1d'
    return "\nError of fixed points array"
  nbins = int(nbins)
  max_f = np.amax(fps)
  min_f = np.amin(fps)
  if mm.size == 2 and mm[1]>mm[0] and (max_f-min_f)>1.0e-12:
    max_f = mm[1]
    min_f = mm[0]
  elif mm.size==1 or mm.size>2:  
    print "\nError: Wrong mm: max and min array!!!"
    return 0
  width = (max_f - min_f)/(2.*nbins)
  if width <1.0e-12:
    hist = np.zeros((2,2))
    hist[0, 0] = max_f
    hist[1, 0] = max_f
    hist[1, 1] = 1.0
    return hist
  hist = np.zeros((nbins, 2))
  hist[:, 0] = width + width*2.0*np.arange(nbins)
  for i in range(fps.size):
    n = int((fps[i] - min_f)/2.0/width)
    if n == nbins:
      n = -1
    elif n>nbins:
      print '\nError: mm[1] is too small!!!'
      print 'mm[1]=', mm[1], ';\nbut max(fps)=', np.amax(fps), '\n'
      return 
    elif n<0:
      print '\nError: mm[0] is too big!!!'
      print 'mm[0]=', mm[1], ';\nbut min(fps)=', np.amin(fps), '\n'
      return 
    hist[n, 1] = hist[n, 1] + 1.0
  if density:
    hist[:,1] = hist[:,1]/1.0/fps.size
  index = find_last_nonzero_hist(hist)
  return hist[:index, :]

# Find the histogram of the array of Fixed Points in a tanh(k) scale
# The histogram is always on the scale of (0, 1)
def hist_tanh_fx(fps, nbins = 100, density = True):
  if fps.ndim > 1:
    print 'Array of fixed point has to be 1d'
    return "\nError of fixed points array"
  nbins = int(nbins)
  fps = np.tanh(fps)
  max_f = np.amax(fps)
  min_f = np.amin(fps)
  width = (max_f - min_f)/(2.*nbins)
  if width <1.0e-15:
    hist = np.zeros((2,2))
    hist[0, 0] = max_f
    hist[1, 0] = max_f
    hist[1, 1] = 1.0
    return hist
  else:
    max_f = 1
    min_f = 0
    width = (max_f - min_f)/(2.*nbins)
    
  hist = np.zeros((nbins, 2))
  hist[:, 0] = width + width*2.0*np.arange(nbins)
  for i in range(fps.size):
    n = int((fps[i] - min_f)/2.0/width)
    if n == nbins:
      n = -1
    elif n>nbins:
      print '\nError: mm[1] is too small!!!'
      print 'mm[1]=', mm[1], ';\nbut max(fps)=', np.amax(fps), '\n'
      return 
    elif n<0:
      print '\nError: mm[0] is too big!!!'
      print 'mm[0]=', mm[1], ';\nbut min(fps)=', np.amin(fps), '\n'
      return 
    hist[n, 1] = hist[n, 1] + 1.0
  if density:
    hist[:,1] = hist[:,1]/1.0/fps.size
  index_1 = find_first_nonzero_hist(hist)
  index_2 = find_last_nonzero_hist(hist)
  return hist[index_1:index_2, :]

# find the last nonzero bins
def find_last_nonzero_hist(hist):
  index = hist.shape[0]-1
  for i in range(hist.shape[0]):
    if hist[index, 1] != 0:
      return index+1
    index = index - 1
  return index

# find the first nonzero bins
def find_first_nonzero_hist(hist):
  for i in range(hist.shape[0]):
    if hist[i, 1] != 0:
      return i
  return i
  
