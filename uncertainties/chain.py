"""

Calculations with full error propagations for quantities described by 
a sample from the underlying distribution. 

Such a sample is a draw from a distribution such that every sample 
is equally likely -- compared to the others -- to be the true value 
of the uncertain quantity (e.g. a MCMC chain).

This module makes several computations based on such samples, 
called chains here, easier to read.

Computations
--------------------------------

For computations with full error propagation with one chain, just
compute values of each sample.

Example:
   
   >>> chain1 = numpy.random.normal(size=1000)
   >>> x = 14 * sin(chain1 * 3 + 2)
   
   The result, x, will be correctly distributed.


"""

import numpy
import scipy.stats
from uncertainties import ufloat, Variable

def mean(chain, **kwargs):
	"""
Obtaining statistical quantities:
-----------------------------------
From the chain, we are interested in things like

  - the mean: chain.mean(chain1) -- numpy.mean
  - the std: chain.std(chain1) -- numpy.std
  - the median: chain.median(chain1) -- numpy.median
  - the quartiles/quantiles: chain.quantiles(chain1, prob=[0.25, 0.75]) -- scipy.stats.mstats.mquantiles
	"""
	return numpy.mean(chain, **kwargs)

def std(chain, **kwargs):
	return numpy.std(chain, **kwargs)
std.__doc__ = mean.__doc__

def median(chain, **kwargs):
	return numpy.median(chain, **kwargs)
median.__doc__ = mean.__doc__

def quantiles(chain, prob, alphap=0, betap=1, **kwargs):
	return scipy.stats.mstats.mquantiles(chain, prob=prob, 
		alphap=alphap, betap=betap, **kwargs)
quantiles.__doc__ = mean.__doc__

def cross(chain1, chain2):
	"""
Combining two independent chains:
-----------------------------------

If chain1 and chain2 were correlated, e.g. one is just a 
linear product of the other, we would do this:

   >>> x = 14 * sin(chain1 * 3 + 2) + 7 * cos(chain2 * 3 + 2)
   
If they are uncorrelated however, we can combine them to the full cross-product.

   >>> C1, C2 = np.meshgrid(chain1, chain2)
   >>> X = 14 * sin(C1 * 3 + 2) + 7 * cos(C2 * 3 + 2)
   >>> x = X.flatten()
   
   or with the syntax from this module:
   
   >>> C1, C2 = chain.cross(chain1, chain2)
   >>> x = chain.combine( 14 * sin(C1 * 3 + 2) + 7 * cos(C2 * 3 + 2) )

This has the benefit that no accidental correlations are introduced, and the 
resulting chain is larger.
"""
	return numpy.meshgrid(chain1, chain2)

def combine(crosschain):
	return crosschain.flatten()

combine.__doc__ = cross.__doc__

def collapse(chain):
	"""
Obtaining a uncertainty:
-----------------------------------
Collapses the chain to a gaussian description (a ufloat).
Be warned of the information loss and that a gaussian may be a poor description
for the distribution at hand.

  >>> v = chain.collapse(chain1)
  >>> print v
"""
	return Variable(mean(chain), std(chain))

class Chain(Variable):
	def __init__(self, chain):
		self.chain = chain
		self.std = std(self.chain)
		self.mean = mean(self.chain)
		super(Chain, self).__init__(self.mean, self.std)
		self.percentiles = quantiles(self.chain, prob=range(101))
	
	def inverse_cdf(self):
		"""
		Returns a function that takes a value between 0 and 1 and 
		returns the corresponding x-value of the cumulative distribution.
		"""
		a = list(self.chain)
		a.sort()
		a = numpy.array(a)
		yinc = a.cumsum()
		x = numpy.linspace(0, 1, len(a))
		return scipy.interpolate.interp1d(x=[0] + list(yinc) + [1], y=[x[0]] + list(x) + [x[-1]])
	
	def __str__(self):
		i = int(numpy.floor(-numpy.log10(self.std_dev()) + 1))
		
		v, errmin, errplus = (self.percentiles[50], self.percentiles[10], self.percentiles[90])
		if i > 0:
			fmt = "%%.%df" % (i)
			return "%s (+%s,-%s)" % (fmt % v,fmt % errplus, fmt % errmin)
		else: # integer representation
			return "%d (+%d,-%d)" % (numpy.round(v, i),numpy.round(errplus, i), numpy.round(errmin, i))
		
	def tex(self):
		i = int(numpy.floor(-numpy.log10(self.std_dev()) + 1))
		
		v, errmin, errplus = (self.percentiles[50], self.percentiles[10], self.percentiles[90])
		if i > 0:
			fmt = "%%.%df" % (i)
			return "{%s}^{+%s}_{-%s}" % (fmt % v,fmt % errplus, fmt % errmin)
		else: # integer representation
			return "{%d}^{+%d}_{-%d}" % (numpy.round(v, i),numpy.round(errplus, i), numpy.round(errmin, i))


__all__ = [collapse, combine, cross, median, std, mean, quantiles]
__doc__ += "\n".join([median.__doc__, cross.__doc__, collapse.__doc__])


