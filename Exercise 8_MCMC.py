# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 10:04:17 2015

@author: Glaucia
"""

"""
Exercise 8 - Introduction to Markov chain Monte Carlo (MCMC)
@author: jembrown
In this example, we will set up a basic MCMC sampler to estimate the probability of success (p) for a binomial distribution.
"""

testPrior = False

from scipy.stats import binom
from scipy.stats import uniform
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

# Define data set
n=20
# increasing the sample size
n = 2000
#decreasing the sample size
n=4
###############################################################################
p_true = 0.5                # True value to be estimated later
data = binom.rvs(n,p_true)  # Drawing a 'data set'
# Define parameter to be estimated and starting value
p = 0.1     # Starting far away from true value to demonstrate burn-in
#Ideally you want to start at different points in sample space to not get stuck in one false top
# Define prior on parameters of interest
def prior(p_in):
    """
    A function to return the prior density for a given value of p.
    """
    return uniform.pdf(p_in,0,1)#giving me the pdf of any point in space, the area under the curve needs to integrate to 1
#The prior is just mantaining us between 0 and 1.
# Define likelihood
def likelihood(p_in):
    """
    A function to return the likelihood for a given value of p.
    """
    if (testPrior):     # Allows us to test the prior, by ignoring the data #if testPrior is true
        return 1
    elif (not testPrior):   # Returns a likelihood score for actual data #if testPrior is false
        if (p_in >= 0 and p_in <= 1):
            return binom.pmf(data,n,p_in)#calculating the laikelihood for different values of p
        else:
            return 0    # Returns value of 0 (rather than non) if p outside bounds.
    
# Define proposal distribution
def drawP(p_in):
    """
    This function provides proposed values for p, given a current value.
    """
    # here a changed the size of the proposal window, by changing the scale values
    newP = uniform.rvs(loc=p_in-0.25,scale=0.5)#scale=10#scale=0.001) # min = loc, max = loc+scale
    return newP
 #In the window of proposed values of p, it is good to mantain the simetry.
#We we want to set a minimum value and the scale is the width of the window.
# Set chain length and run it
ngens = 50          # Total length of chain
sampleFreq = 1      # Change this to something >1 if you want to space out samples.
updateFreq = 5000   #how often this should print results on the screen # Frequency of screen updates to make sure chain is running.
samples = []        # Vector to hold sampled values
# recording samples for each generation:
ngens = 50          # Total length of chain
sampleFreq = 1      # Change this to something >1 if you want to space out samples.
updateFreq = 50     # Frequency of screen updates to make sure chain is running.
samples = []        # Vector to hold sampled values
#spacing out samples
ngens = 50          # Total length of chain
sampleFreq = 10      # Change this to something >1 if you want to space out samples.
updateFreq = 5000   # Frequency of screen updates to make sure chain is running.
samples = []        # Vector to hold sampled values

for gen in range(ngens):
    p_prop = drawP(p)
    p_prop_post = prior(p_prop)*likelihood(p_prop)
    p_post = prior(p)*likelihood(p)
    if (p_prop_post >= p_post): # If proposed value has posterior density > curr value
        p = p_prop
    elif (p_prop_post < p_post): # If proposed value has posterior density < curr value
        r = p_prop_post/p_post#how most worse it is in relation to your current value
        ranUnifDraw = uniform.rvs()#draw a random number from an unifrom and compare to r
        if (ranUnifDraw <= r):
            p = p_prop
    else:
        print "Problem calculating proposal ratio!"
        print (p_prop,p_prop_post)
        print (p,p_post)
        print ""
    if (gen % sampleFreq == 0):#is doing division, but instead of giving the result, it gives the remainder
        samples.append(p)
    if (gen % updateFreq == 0):
        print "Generation %s" % gen

# Summarizing MCMC samples

# Numerical summaries
burnin = int((ngens/sampleFreq)*0.1)
postBurnSamples = samples[burnin+1:]#you want to look to your data before defining it. There is a stationary distribution, the posterior is the stationary
print("Posterior Mean: %f" % np.mean(postBurnSamples))#mean of the posterior 
postBurnSamples.sort()  # post-burnin samples will be sorted after this is called # you need to sort it, otherwise it will give you wrong values
print("Posterior 95% credibility interval: "+"(%f,%f)" % (postBurnSamples[int(len(postBurnSamples)*0.025)],postBurnSamples[int(len(postBurnSamples)*0.975)]))
# In Bayesian - Credibel Intervals (not Confidence Intervals)
# Marginal histogram
pl.figure()
pl.hist(postBurnSamples)
pl.show()

# Trace plot
plt.figure()
plt.plot(range(ngens/sampleFreq),samples)
plt.ylim(0,1)
plt.ylabel("P")
plt.xlabel("Sample Number")
plt.show()      

"""
Play around with different settings for the analyses outlined above and try to 
answer these questions:
(1) How do trace plots differ if we record samples from every generation or if 
we space out our sampling?
"""
##If we record samples from every generation or if we space out the sampling, the trace plot will look less precise. For example, if the trace was a wave, the wave would have lower frequency, with bigger steps between valleys and peaks.. 
"""
(2) How does the width of the 95% credible interval change as we add/subtract data?
"""
##As we add data, the widht of th 95% credible interval decreases. As we subtract data, its width incresases.
"""
(3) How does the trace plot change if the proposal window becomes much bigger 
or smaller for a given amount of data?
"""
##With a small amount of data the trace plot gets more stable around the posterior mean, while with small ammount of data, the trace plot is more unstable and differs more often from the mean.
##With the same ammount of data, but varying the proposal window. Much bigger: the trace plot becomes constant in some points, but can deviate abruptly from the mean.
##Much smaller: trending appear on the trace plot