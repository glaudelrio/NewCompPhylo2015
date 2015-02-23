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
def BinomSimp(n,k):##Pated this function to calculate the factorial in the BinomPMF function below.
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
   
def binomPMF(n,k,p):##This function calculates the binomial PMF
    binomcoef=BinomSimp(n,k)##factorial of n and k
    return binomcoef*(pow(p,k))*(pow(1-p,n-k))##multiplying the factorial by the rest of the equation of binomial probability mass function

binomPMF(5,4,0.33)##testing the function

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
print(pvalues)##a list with values of p going from 0 to 1 by steps of 0.05
# Calculate the likelihood scores for these values of p, in light of the data you've collected
Lik=[]
for prob in pvalues:
    liks=binomPMF(5,4,p=prob)##as I have the number of k, I will calculate the likelihhod of each value of p based on the number os successes obtained (k), I can do this with the binomial PMF
    Lik.append(liks)##putting the likelihhods in a list
print(Lik)
# Find the maximum likelihood value of p (at least, the max in this set)
print max(Lik)##max likelihood associated to p=0.8
# What is the strength of evidence against the most extreme values of p (0 and 1)?
#pA(x)/pb(x)
0.4096/0## as the Likelihoods for p=0 or p=1 are 0, the strength of evidence is undefined
# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
LikRatios=[]
for x in Lik:
    division=x/0.4096##to calculate the LRatios, I just divided each likelihood (for each p) by the maximum likelihhod
    LikRatios.append(division)##putting the LRatios in a list
print(LikRatios)
#################Trying to win this bet:
values2=range(70,90,1)##here I changed the intervals to see if I would got the same maximum likelihhod estimator as before
pvalues2=[]
for i in values2:
    result2=i/100
    pvalues2.append(result2)
print(pvalues2)
Lik2=[]
for prob in pvalues2:
    liks2=binomPMF(5,4,p=prob)
    Lik2.append(liks2)
print(Lik2)
max(Lik2)
#################And one more time
values3=range(795,805,1)##here I used an interval even more precise
pvalues3=[]
for i in values3:
    result3=i/1000
    pvalues3.append(result3)
print(pvalues3)
Lik3=[]
for prob in pvalues3:
    liks3=binomPMF(5,4,p=prob)
    Lik3.append(liks3)
print(Lik3)
max(Lik3)####associated with 0.8
"""
Now let's try this all again, but with more data. This time, we'll use 20 draws from our cup of marbles.
"""
data = 12 # Supply observed number of successes here.
numTrials = 20
# Calculate the likelihood scores for these values of p, in light of the data you've collected
values4=range(20,80,5)##used an interval of probs going from 0.2 to 0.8 
pvalues4=[]
for i in values4:
    result4=i/100
    pvalues4.append(result4)
print(pvalues4)
Lik4=[]
for prob in pvalues4:
    liks4=binomPMF(20,12,p=prob)
    Lik4.append(liks4)##list of likelihhods associated to each value of p for k=12
print(Lik4)
# Find the maximum likelihood value of p (at least, the max in this set)
max(Lik4)###associated with p=0.6
# What is the strength of evidence against the most extreme values of p (0 and 1)?
####Again the force of evidence is undefined, since the likelihood for 0 and for 1 is equal to 0.
# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
LikRatios2=[]
for x in Lik4:
    division2=x/0.17970578775468937##Again, calculated the LRatios by dividing each Likelihood by the maximum likelihood
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
def mle(n,k,pCurr,diff):##this function will give you the maximum likelihood estimator, if you give the number of trials (n), number of successes (k),an initial value for p (pCurr) and a diff. The diff represents the size of the steps the function will take when searching the mle. So if you use a small diff (0.0001), you will get a more precise estimate.
    pUp=pCurr+diff
    pDown=pCurr-diff
    vpCurr=binomPMF(n,k,p=pCurr)##Likelihhod of current p
    vpUp=binomPMF(n,k,p=pUp)##likelihood of pUp (p+diff)
    vpDown=binomPMF(n,k,p=pDown)##likelihhod of pDown (p-diff)
    if (vpCurr < vpUp): ##if the likeliood of vpCurr is small than the likelihood of vpUp go to the next loop   
        while (vpCurr < vpUp):
            pCurr=pUp##this will make the function "walk" to the direction of the mle
            vpCurr=binomPMF(n,k,p=pCurr)##calculate the likelihood of current p
            pUp=pCurr+diff##add diff to the current p
            vpUp=binomPMF(n,k,p=pUp)##calculate the likelihood for current p + diff ##while this values is bigger than the value associated to the current p, the loop continues
        return pCurr##when the maximum likelihhod is reached, it will return the mle
    else:##same as above but to the opposite direction
        while (vpCurr < vpDown):
            pCurr=pDown
            vpCurr=binomPMF(n,k,p=pCurr)
            pDown=pCurr-diff
            vpDown=binomPMF(n,k,p=pDown)
        return pCurr
    
##testing the function,
banana=mle(5,4,0.9,0.0001)
        
print(banana)
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
for i in range (10):
    outcome=DiscreteSample(xk=[0,1],pk=[0.3,0.7],siz=200)
    counting1=outcome.count(1)
    sim.append(counting1)
print sim
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
print LtrueP

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
        return vpUp
    else:
        while (vpCurr < vpDown):
            pCurr=pDown
            vpCurr=binomPMF(n,k,p=pCurr)
            pDown=pCurr-diff
            vpDown=binomPMF(n,k,p=pDown)
        return vpDown

ML=[]
for y in sim:
    estimate=ml(200,k=y,pCurr=0.01,diff=0.001)
    ML.append(estimate)
print(ML)
#########Getting the Ratios
LRatios = [float(b) / float(m) for b,m in zip(LtrueP, ML)]
print (LRatios)
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
# Any idea why that particular cell would be meaningful?
# Based on your results (and the values in the table), what LR statistic value
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?
# Using this cutoff, what interval might you report for the 5- and 20-trial data
# sets above?
# We've talked in previous classes about two ways to interpret probabilities. Which
# interpretation are we using here to define these intervals?