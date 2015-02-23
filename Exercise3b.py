

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 01 18:01:25 2015
@author: Glaucia
"""

"""
An Introduction to Likelihood
@author: jembrown
"""
"""
There are two primary ways to use probability models. Given what we know to be true about an experimental setup, we can make predictions about what we expect to see in an upcoming trial. For this purpose, probability functions are what we need. If R is the outcome of an experiment (i.e., an event) and p is a parameter of the probability function defined for these experimental outcomes, we can write these expectations or predictions as:
P(R|p).
This conditional probability tells us what to expect about the outcomes of our experiment, given knowledge of the underlying probability model. Thus, it is a function of R given a value for p (and the model itself).
However, we might also wish to ask what we can learn about p, given outcomes of trials that have already been observed. This is the purview of the likelihood. Likelihoods are functions of parameters (or hypotheses, more generally) given some observations. The likelihood function of a parameter value is defined as:
L(p;R) = P(R|p)
This is the same probability statement we saw above. However, in this context we are considering the outcome (R) to be fixed and we're interested in learning about p. Note that the likelihood is sometimes written in several different ways: L(p;R) or L(p) or L(p|R). P(R|p) gives a probability when R is discrete or a probability density when R is continuous. Since likelihoods are only compared for some particular R, we do not need to worry about this distinction. Technically speaking, likelihoods are just said to be proportional to P(R|p), with the constant of proportionality being arbitrary.
There are some very important distinctions between likelihoods and probabilities. First, likelihoods do NOT sum (or integrate) to 1 over all possible values of p. Therefore, the area under a likelihood curve is not meaningful, as it is for probability.
It does not make sense to compare likelihoods across different R. For instance, smaller numbers of observations generally produce higher values of P(R|p), because there are fewer total outcomes.
Likelihood curves provide useful information about different possible values of p. When we are interested in comparing different hypotheses (H1 and H2), the likelihood ratio is often used:
L(H1;R) P(R|H1)
------- = -------
L(H2;R) P(R|H2)
Now, let's try using likelihoods to learn about unknown aspects of the process that's producing some data.
---> Inferring p for a binomial distribution <---
First, we'll start by trying to figure out the unknown probability of success associated with a Binom(5,p) random variable. If you want to try this on your own later, the following code will perform draws from a binomial with 5 trials. You can simply change the associated value of p to whatever you'd like. To make the inference blind, have a friend set this value and perform the draws from the Binomial for you, without revealing the value of p that they used.
"""
"""
For the in-class version of this exercise, I'm going to perform a manual draw from a binomial using colored marbles in a cup. We'll arbitrarily define dark marbles as successes and light marbles as failures.
Record the outcomes here:
Draw 1:Dark
Draw 2:Dark
Draw 3:Dark
Draw 4:Dark
Draw 5:White
Number of 'successes':
Now record the observed number of succeses as in the data variable below.
"""
data = 4 # Supply observed number of successes here.
numTrials = 5
"""
Since we are trying to learn about p, we define the likelihood function as;
L(p;data) = P(data|p)
If data is a binomially distributed random variable [data ~ Binom(5,p)]
P(data=k|p) = (5 choose k) * p^k * (1-p)^(n-k)
So, we need a function to calculate the binomial PMF. Luckily, you should have just written one and posted it to GitHub for your last exercise. Copy and paste your binomial PMF code below. For now, I will refer to this function as binomPMF().
"""
def BinomSimp(n,k):
   nb1=1
   kb1=1
   sub = n-k
   while (n > sub):
       nb1 = nb1*n
       n = n-1
   while (k > 0):
       kb1 = kb1*k
       k= k-1
       result=nb1/kb1
   return result
   
def binomPMF(n,k,p):
    binomcoef=BinomSimp(n,k)
    return binomcoef*(pow(p,k))*(pow(1-p,n-k))

binomPMF(5,4,0.33)

"""
Now we need to calculate likelihoods for a series of different values for p to compare likelihoods. There are an infinite number of possible values for p, so let's confine ourselves to steps of 0.05 between 0 and 1.
"""
# Set up a list with all relevant values of p
from __future__ import division
values=range(0,105,5)
pvalues=[]
for i in values:
    result=i/100
    pvalues.append(result)
print(pvalues)
# Calculate the likelihood scores for these values of p, in light of the data you've collected
Lik=[]
for prob in pvalues:
    liks=binomPMF(5,4,p=prob)
    Lik.append(liks)
print(Lik)
# Find the maximum likelihood value of p (at least, the max in this set)
print max(Lik)#max likelihood associated to p=0.8
# What is the strength of evidence against the most extreme values of p (0 and 1)?
#pA(x)/pb(x)
0.4096/0# as the Likelihoods for p=0 or p=1 are 0, the strength of evidence is undefined
# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
LikRatios=[]
for x in Lik:
    division=x/0.4096
    LikRatios.append(division)
print(LikRatios)

"""
Now let's try this all again, but with more data. This time, we'll use 20 draws from our cup of marbles.
"""
data = 12 # Supply observed number of successes here.
numTrials = 20
# Calculate the likelihood scores for these values of p, in light of the data you've collected
values4=range(20,80,5)
pvalues4=[]
for i in values4:
    result4=i/100
    pvalues4.append(result4)
print(pvalues4)
Lik4=[]
for prob in pvalues:
    liks4=binomPMF(20,12,p=prob)
    Lik4.append(liks4)
print(Lik4)
# Find the maximum likelihood value of p (at least, the max in this set)
max(Lik4)###associated with p=0.6
# What is the strength of evidence against the most extreme values of p (0 and 1)?
####Again the force of evidence is undefined, since the likelihood for 0 and for 1 is equal to 0.
# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
LikRatios2=[]
for x in Lik4:
    division2=x/0.17970578775468937
    LikRatios2.append(division2)
print(LikRatios2)
# When is the ratio small enough to reject some values of p?
####You could consider that the ratio is small enough when it is less than 0.05 (Savage et al. 2011)
# Note: You will empirically investigate this on your own later in this exercise.
# **** EVERYTHING ABOVE HERE TO BE POSTED TO GITHUB BY TUESDAY, FEB. 3RD. ****
# **** CODE BELOW TO BE POSTED TO GITHUB BY THURSDAY, FEB. 5TH ****
"""
Sometimes it will not be feasible or efficient to calculate the likelihoods for every
value of a parameter in which we're interested. Also, that approach can lead to large
gaps between relevant values of the parameter. Instead, we'd like to have a 'hill
climbing' function that starts with some arbitrary value of the parameter and finds
values with progressively better likelihood scores. This is an ML optimization
function. There has been a lot of work on the best way to do this. We're going to try
a fairly simple approach that should still work pretty well, as long as our likelihood
surface is unimodal (has just one peak). Our algorithm will be:
(1) Calculate the likelihood for our starting parameter value (we'll call this pCurr)
(2) Calculate likelihoods for the two parameter values above (pUp) and below (pDown)
our current value by some amount (diff). So, pUp=pCurr+diff and pDown=pCurr-diff. To
start, set diff=0.1, although it would be nice to allow this initial value to be set
as an argument of our optimization function.
(3) If either pUp or pDown has a better likelihood than pCurr, change pCurr to this
value. Then repeat (1)-(3) until pCurr has a higher likelihood than both pUp and
pDown.
(4) Once L(pCurr) > L(pUp) and L(pCurr) > L(pDown), reduce diff by 1/2. Then repeat
(1)-(3).
(5) Repeat (1)-(4) until diff is less than some threshold (say, 0.001).
(6) Return the final optimized parameter value.
Write a function that takes some starting p value and observed data (k,n) for a
binomial as its arguments and returns the ML value for p.
To write this function, you will probably want to use while loops. The structure of these loops is
while (someCondition):
code line 1 inside loop
code line 2 inside loop
As long as the condition remains True, the loop will continue executing. If the
condition isn't met (someCondition=False) when the loop is first encountered, the
code inside will never execute.
If you understand recursion, you can use it to save some lines in this code, but it's
not necessary to create a working function.
"""
# Write a function that finds the ML value of p for a binomial, given k and n.
def mle(n,k,pCurr,diff):##
    """
    I tested bigger and smaller values for diff, as you mentioned in the statements above and decided to live diff as an argument. If you use small values, you will get a better estimate.
    """
    pUp=pCurr+diff
    pDown=pCurr-diff
    vpCurr=binomPMF(n,k,p=pCurr)
    vpUp=binomPMF(n,k,p=pUp)
    vpDown=binomPMF(n,k,p=pDown)
    if (vpCurr < vpUp):    
        while (vpCurr < vpUp):
            pCurr=pUp
            vpCurr=binomPMF(n,k,p=pCurr)
            pUp=pCurr+diff
            vpUp=binomPMF(n,k,p=pUp)
        return pCurr
    else:
        while (vpCurr < vpDown):
            pCurr=pDown
            vpCurr=binomPMF(n,k,p=pCurr)
            pDown=pCurr-diff
            vpDown=binomPMF(n,k,p=pDown)
        return pCurr
    

mle(5,4,0.9,0.0001)
        

"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of
p. Now, we will empirically determine one way to construct such an interval. To do
so, we will ask how far away from the true value of a parameter the ML estimate
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once
you have this distribution, find the likelihood ratio cutoff you need to ensure
that the probability of seeing an LR score that big or greater is <= 5%.
"""
# Set a starting, true value for p
trueP =0.70
# Simulate 1,000 datasets of 200 trials from a binomial with this p
from scipy.stats import rv_discrete

def DiscreteSample (xk,pk,siz):
    discrete = rv_discrete(name='discrete', values=(xk, pk))
    sample = discrete.rvs(size=siz)
    x = list(sample)
    return x

sim=[]
for i in range (1001):##simulating 1001 datasets
    outcome=DiscreteSample(xk=[0,1],pk=[0.3,0.7],siz=200)##with 200 trials
    counting1=outcome.count(1)##counting the frequencies of one of the outcomes
    sim.append(counting1)##appending the frequencies in a list
print sim##printing the list
# If you haven't already done so, you'll want to import the binom class from scipy:
# from scipy.stats import binom
# binom.rvs(n,p) will then produce a draw from the corresponding binomial.
# Now find ML parameter estimates for each of these trials
MLE=[]
for y in sim:
    estimate=mle(200,k=y,pCurr=0.01,diff=0.001)
    MLE.append(estimate)
print(MLE)

#to be continue...
# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum
# likelihood (ML) in the denominator. Sort the results and find the value
# corresponding to the 95th percentile.
####Getting the L(trueP):
LtrueP=[]
for t in sim:
    Ltp=binomPMF(200,k=t,p=0.7)
    LtrueP.append(Ltp)
print LtrueP####List with Ltruep values

######Getting the ML values:
def ml(n,k,pCurr,diff):
    pUp=pCurr+diff
    pDown=pCurr-diff
    vpCurr=binomPMF(n,k,p=pCurr)
    vpUp=binomPMF(n,k,p=pUp)
    vpDown=binomPMF(n,k,p=pDown)
    if (vpCurr < vpUp):    
        while (vpCurr < vpUp):
            pCurr=pUp
            vpCurr=binomPMF(n,k,p=pCurr)
            pUp=pCurr+diff
            vpUp=binomPMF(n,k,p=pUp)
        return vpCurr
    else:
        while (vpCurr < vpDown):
            pCurr=pDown
            vpCurr=binomPMF(n,k,p=pCurr)
            pDown=pCurr-diff
            vpDown=binomPMF(n,k,p=pDown)
        return vpCurr

ml(5,4,0.2,diff=0.01)#testing the function

ML=[]
for y in sim:
    estimate=ml(200,k=y,pCurr=0.01,diff=0.01)
    ML.append(estimate)
print(ML)#####List with ML values
#########Getting the Ratios
LRatios = [float(b) / float(m) for b,m in zip(LtrueP, ML)]
print (LRatios)
#########Getting the number corresponding to the 95th percentile
from numpy import percentile
print percentile(LRatios,95)

# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values.
from math import log
neg2ln=[]
for w in LRatios:
    loglikneg=(-2)*(log(w))
    neg2ln.append(loglikneg)
print neg2ln
    
# Find the 95th percentile of these values. Compare these values to this table:# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look
# at the 0.05 column. Do any of these values seem similar to the one you calculated?
print percentile(neg2ln,95)
#########3.85
#########Yeah, sure! It is very similar to the number in the cell representing 1 degree of freedom and alfa of 0.05...
###############################################################################
# Any idea why that particular cell would be meaningful?
##########This cell could be meaningful because in a binomial trial we are dealing with only one degree of freedom...
###############################################################################
# Based on your results (and the values in the table), what LR statistic value
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?
#########It seems that values above that 95% percentile of [-2ln(LR)] will indicate the less probable values of p (appearing only in 5% of the cases)
###############################################################################
# Using this cutoff, what interval might you report for the 5- and 20-trial data
# sets above?
#########5 trials
#########For this I will have to calculate the -2ln(LR)

LikRatios.pop(0)###I had to take this value from my function to make the log work
LikRatios.pop(19)###I had to take this value from my function to make the log work
stat1=[]
for b in LikRatios:
    banana=(-2)*(log(b))
    stat1.append(banana)
print stat1
#########Using 3.85 as cutoff I was able to exclude p=0 to p=0.35 and also p=1. So now I have a CI in between 0.40 to 0.95.I did it visually, just scanning my list.
#########20 trials
LikRatios2.pop(0)
LikRatios2.pop(19)
stat2=[]
for t in LikRatios2:
    picles=(-2)*(log(t))
    stat2.append(picles)
print stat2
#########Using 3.85 as cutoff I was able to exclude p=0 to p=0.35 and also p=0.8 to p=1. So now I have a CI in between 0.40 to 0.75.
###############################################################################
# We've talked in previous classes about two ways to interpret probabilities. 
# Which interpretation are we using here to define these intervals?
#########This interpretation is based on a test of frequencies. It could be seen as a frequentist approach.
###############################################################################