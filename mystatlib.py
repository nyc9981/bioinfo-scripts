# Could use iterable as the name of parameter instead of arr
import math
import numpy as np
from collections import Counter

def mean(arr):
    return sum(arr) / float(len(arr))

def mean_recursive(arr):
    """ Recursive implementation of mean function
        Prevents overflow
    """
    n = len(arr)
    if n == 1:
        return arr[0]
    return mean_recursive(arr[:n-1]) * (1-1./n) + arr[-1]/float(n)

def variance(arr):
    """ Population Variance
    """
    X2 = [x**2 for x in arr]
    return mean(X2) - mean(arr)**2

def std_dev(arr): 
    """ Population Standard Deviation
    """
    return variance(arr)**0.5

def var(arr):
    """ Sample Variance
    """
    return variance(arr) / (1 - 1.0 / len(arr))

def sd(arr):
    """ Sample Standard Deviation
    """
    return var(arr)**0.5

def median(arr):
    sorted_arr = sorted(arr)
    n = len(sorted_arr) 
    return sorted_arr[n//2] if n%2 != 0 else (sorted_arr[n//2] + sorted_arr[n//2 -1])/2.

# from collections import Counter
def mode(arr):
    cer = Counter(arr).most_common()
    return sorted(cer, key=lambda (k, v): (-v, k))[0][0] # mode (minimum to break ties)

def quartiles(arr):
    """ Calculate the 1st, 2nd, 3rd Quartiles
    """
    def median(sorted_arr):
        n = len(sorted_arr)
        return sorted_arr[n//2] if n%2 != 0 else (sorted_arr[n//2] + sorted_arr[n//2 -1])/2.
    n = len(arr)
    a = sorted(arr)
    q2 = median(a)
    L = a[: n//2] # lower half
    U = a[n//2 :] if n%2==0 else a[n//2+1 :] # upper half
    q1 = median(L)
    q3 = median(U)
    return (q1, q2, q3)

def iqr(ele, frq):
    """ Interquartile Range of data that contain frequencies of each element
    """
    n = len(ele)
    arr = sorted(zip(ele, frq))
    sorted_ele = [arr[i][0] for i in xrange(n)]
    sorted_frq = [arr[i][1] for i in xrange(n)]
    
    sorted_arr = []
    for i in xrange(n):
        sorted_arr += [sorted_ele[i]]*sorted_frq[i]
    q1, _, q3 = quartiles(sorted_arr)
    return q3 - q1

def cov(X, Y):
    """ Covariance
    """ 
    nx, ny = len(X), len(Y)
    if nx != ny:
        #return None # raise some knid of exception?
        raise ValueError('X and Y must have the same length')
   
    t = sum((X[i] - X[j]) * (Y[i] - Y[j]) for i in xrange(nx) for j in xrange(i+1, nx) )
    return 1. / nx**2 * t  

def pearson_corr_coef(X, Y):
    """ Person's Correlation Coefficienct
    """ 
    return cov(X, Y) / std_dev(X) / std_dev(Y)

def rank(X):
    """ Rank in ascending order
    """ 
    uniq_sorted_X = sorted(list(set(X)))
    return [uniq_sorted_X.index(x) + 1 for x in X]

def spearman_rank_corr_coef(X, Y):
    """ Spearman's Rank Correlation Coefficienct
    """ 
    X_ranks = rank(X)
    Y_ranks = rank(Y)
    return pearson_corr_coef(X_ranks, Y_ranks)

class RunningStat:
    '''
    The class RunningStat uses its methods to compute the mean, 
    sample variance, and standard deviation of a stream of data.
    See https://www.johndcook.com/blog/standard_deviation/
    '''
    def __init__(self, m_n = 0):
        self.m_n = m_n
        
    @classmethod    
    def from_array(cls, arr):
        rs = cls()
        for a in arr:
            rs.push(a)
        return rs
    
    def clear(self):
        self.m_n = 0
        
    def push(self, x):
        self.m_n += 1
        
        # See Knuth TAOCP vol 2, 3rd edition, page 232
        if self.m_n == 1:
            self.m_oldM = self.m_newM = x
            self.m_oldS = 0.0
        else:
            self.m_newM = self.m_oldM + (x - self.m_oldM)/float(self.m_n)
            self.m_newS = self.m_oldS + (x - self.m_oldM)*(x - self.m_newM)
    
            # set up for next iteration
            self.m_oldM = self.m_newM; 
            self.m_oldS = self.m_newS;
 
    def num_data_values(self):
        return self.m_n

    def mean(self):
        return self.m_newM if self.m_n > 0 else 0.0

    def var(self):
        """ Sample Variance: Estimator of Population Variance
        """
        return self.m_newS/(self.m_n - 1) if self.m_n > 1 else 0.0 

    def sd(self):
        """ Sample Standard Deviation
        """
        return self.var() ** 0.5        
        
    def variance(self):
        """ Population Variance
        """
        return self.m_newS/self.m_n if self.m_n > 0 else 0.0 

    def std_dev(self):
        """ Population Standard Deviation
        """
        return self.variance() ** 0.5

# Simple Linear Regression
def LR(X, Y): 
    """ Linear Regression with one feature in addition to 1's 
        Least Square Regression
        Returns a predict function
    """
    b = pearson_corr_coef(X, Y) * std_dev(Y) / std_dev(X)
    a = mean(Y) - b * mean(X)
    return lambda x: a + b*x

# Multiple Linear Regression
# import numpy as np
def MLR(X, Y):
    """ Linear Regression with multiple features in addition to 1's
        Multiple Least Square Regression
        Returns a predict function
    """
    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    return lambda x: B.dot(x.T) 

#{{{{{{{{---------------- PROBABILITY DISTRIBUTIONS ----------------
def choose(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def dbinom(k, n, p):
    """ Binomial Distribution: Probability Mass Function (PMF)
        k successes, given n trials
    """
    return choose(n, k) * p**k * (1-p)**(n-k)

def pbinom(k, n, p):
    """ Binomial distribution: Cumulative Distribution Function (CDF)
        k successes, given n trials
    """
    return sum(dbinom(x, n, p) for x in xrange(k+1))


def dnbinom(n, k, p):
    """ Negative Binomial Distribution: Probability Mass Function (PMF)
        n trials or r=n-k failures, given k successes
    """
    return choose(n-1, k-1) * p**k * (1-p)**(n-k)

def pnbinom(n, k, p):
    """ Negative binomial distribution: Cumulative Distribution Function (CDF)
        n trials or r=n-k failures, given k successes
    """
    return sum(dnbinom(x, k, p) for x in xrange(1, n+1))


def dgeom(n, p):
    """ Geometric Distribution: Probability Mass Function (PMF)
        n trials, given 1st success
    """
    return dnbinom(n, 1, p)

def pgeom(n, p):
    """ Geometric distribution: Cumulative Distribution Function (CDF)
        n trials, given 1st success
    """
    return pnbinom(n, 1, p)

def dhgeom(k, w, b, n):
    """ Hypergeometric distribution: Probability Mass Function (PMF)
        k successes (desired objects) given n draws from w desired and b undesired objects
    """
    return choose(w, k) * choose(b, n-k) / choose(w+b, n)

def phgeom(k, w, b, n):
    """ Hypergeometric distribution: Cumulative Distribution Function (CDF)
        k successes (desired objects) given n draws from w desired and b undesired objects
    """
    return sum(dhgeom(x, w, b, n) for x in xrange(k+1))


#import math
def dpois(k, ld):
    """ Poisson distribution: Probability Mass Function (PMF)
        ld: {lambda}
    """
    return ld ** k * math.e ** -ld / math.factorial(k) 

def ppois(k, ld):
    """ Poisson distribution: Cumulative Distribution Function (CDF)
        ld: {lambda}
    """
    return sum(dpois(i, ld) for i in xrange(k+1))

#import math
def pnorm(x, u, s):
    """ Normal distribution: Cumulative Distribution Function (CDF)
    """    
    return 1./2 * (1 + math.erf( (x - u) / (s * math.sqrt(2)) ))

#}}}}}}}}---------------- PROBABILITY DISTRIBUTIONS ----------------

def main():
    #print "HELLO"
    a = [1, 3., 5, 2.2, 5, 9, 10.2, 10.1, 7.8, 9.3, 13.2]
    print mean(a), variance(a), std_dev(a), var(a), sd(a)
    
    rs = RunningStat()
    for x in a:
        rs.push(x)
    print rs.mean(), rs.variance(), rs.std_dev(), rs.var(), rs.sd()
    
    rs2 = RunningStat.from_array(a)
    print rs2.mean(), rs2.variance(), rs2.std_dev(), rs2.var(), rs2.sd()
    
    print mode(a)

if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    #test
    main()