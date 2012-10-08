import numpy

from uncertainties import chain

def test_runthrough():
	chain1 = numpy.random.normal(size=100)
	x = 14 * numpy.sin(chain1 / 100 + 2)

	print "mean     ", chain.mean(x)
	print "std      ", chain.std(x)
	print "median   ", chain.median(x)
	print "quantiles", chain.quantiles(x, [0.25, 0.5, 0.75])

	c = chain.Chain(x)
	assert len(c.percentiles) == 101
	print "string   ", c

	print "full repr", repr(c)
	
	assert c.inverse_cdf()(0.5) == c.median, [c.inverse_cdf()(0.5), c.median]

def test_2chains():
	chain1 = numpy.random.normal(size=100)
	chain2 = numpy.random.normal(size=100)

	C1, C2 = chain.cross(chain1, chain2)
	x = chain.combine( 14 * numpy.sin(C1 * 3 + 2) + 7 * numpy.cos(C2 * 3 + 2) )

	print "mean     ", chain.mean(x)
	print "std      ", chain.std(x)
	print "median   ", chain.median(x)
	print "quantiles", chain.quantiles(x, [0.25, 0.5, 0.75])


