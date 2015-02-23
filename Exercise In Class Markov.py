# -*- coding: utf-8 -*-
"""
Glaucia Del-Rio

"""

"""
Recall from your reading that any irreducible and aperiodic Markov chain has a 
stationary distribution. To convince ourselves that things will converge for 
such a chain with arbitrary transition probabilities, let's give it a try.
Work in pairs for this. It's more fun to be social.
"""

# Paste your Markov chain simulation function below, where the starting state
# is drawn with uniform probability from all possible states. Remember to also
# copy any import statements or other functions on which your simulator is
# dependent.
from numpy.random import uniform
mat = [[0.3,0.7],[0.4,0.6]]

# Paste or import your discrete sampling function
import scipy       
def discSamp(events,probs):
    ranNum = scipy.random.random()
    cumulProbs = []
    cumulProbs.extend([probs[0]])
    for i in range(1,len(probs)):
        cumulProbs.extend([probs[i]+cumulProbs[-1]])
    for i in range(0,len(probs)):
        if ranNum < cumulProbs[i]:
            return events[i]
    return None
# Write your Markov chain simulator below. Record the states of your chain in
# a list. Draw a random state to initiate the chain. 
####My states are 0 and 1:
def MarkovChain (x,matrix):
    rn=uniform(low=0.0, high=1.0, size=None)
    currState=discSamp(events=[0,1],probs=[rn,1-rn])
    list1=[currState]    
    for i in range (x):
        currState=discSamp(events=[currState,1-currState],probs=matrix[currState])
        list1.append(currState)
    return list1

MarkovChain (40,matrix=mat)

###############################################################################

####To multiplicate a matrix by other in numpy... you could raise it to a power as:
import numpy
m=numpy.matrix(mat)
m**6
# Define a 2x2 transition matrix. For fun, don't make all the probabilities
# equal. Also, don't use any 0s or 1s (to make sure the chain is irreducible
# and aperiodic).
mat = [[0.3,0.7],[0.4,0.6]]

# Simulate a single chain for three time steps and print the states
MarkovChain (2,matrix=mat)

# Analytically calculate the progression of states for this chain.
####P(A=1)*P(B=0|A=1)*P(C=1|B=0)
Progression=0.5*0.6*0.7
print(Progression)
# Calculate the probability of observing the state in step 3, given the initial
# state in step 1 (i.e., as if you didn't know the state in step 2).
####P(A=1)*P(B=0|A=1)*P(C=1|B=0)+P(A=1)*P(B=1|A=1)*P(C=1|B=1)
Prob1_3=(0.5*0.6*0.7)+(0.5*0.4*0.4)
print(Prob1_3)

# Now think of the chain progressing in the opposite direction. What is the
# probability of the progression through all 3 states in this direction? How
# does this compare to the original direction?
####It has the same value.
RevProgression=0.5*0.6*0.7
print(RevProgression)
####The probability of observing a particular sequence is the same for both directions (Reversibility) 

# Try the same "forward" and "reverse" calculations as above, but with this
# transition matrix:
# revMat = [[0.77,0.23],
#           [0.39,0.61]]
# and these starting frequencies for "a" and "b"
# freq(a) = 0.63   freq(b) = 0.37
####Redefining the function to change the initial frequencies.
def MarkovChain2 (x,matrix):
    currState=discSamp(events=[0,1],probs=[0.63,0.37])
    list1=[currState]    
    for i in range (x):
        currState=discSamp(events=[currState,1-currState],probs=matrix[currState])
        list1.append(currState)
    return list1

mat2=[[0.77,0.23],[0.61,0.39]]

MarkovChain2 (2,matrix=mat2)
####[1,1,0]
####Progression 2
P2=0.37*0.39*0.61
print P2
####RevP2
RevP2=0.63*0.23*0.39
print RevP2

# What is (roughly) true about these probabilities?
####The probabilities change according to the order imposed to the chain, this happened because we set different initial probabilities for the states.

# Simulate 1,000 replicates  (or 10K if your computer is fast enough) of 25 
# steps. What are the frequencies of the 2 states across replicates through time?
list1=[]
for i in range(1000):
    result=MarkovChain2 (25,matrix=mat2)
    list1.append(result)
    
print list1

# NOTE: Here is a function that reports the frequencies of a state through time 
# for replicate simulations. You'll need to do this several times during this exercise.

def mcStateFreqSum(sims,state="a"):
    """
    Pass this function a list of lists. Each individual list should be the
    states of a discrete-state Markov chain through time (and all the same 
    length). It will return a list containing the frequency of one state 
    ("a" by default) across all simulations through time.
    """
    freqs = []
    for i in range(len(sims[0])):  # Iterate across time steps
        stateCount = 0
        for j in range(len(sims)): # Iterate across simulations
            if sims[j][i] == state:
                stateCount += 1
        freqs.extend([float(stateCount)/float(len(sims))])
    return freqs

# Run replicate simulations 

repfreq1=mcStateFreqSum(sims=list1,state=0)

# Summarize the frequency of one state through time
    
print repfreq1

# What do you notice about the state frequencies through time? Try another round
# of simulations with a different transition matrix. How do the state freq.
# values change?
####I noticed that the frequencies converge to a stationary frequency, in this case near by 0.63.
mat3=[[0.23,0.77],[0.39,0.61]]
list2=[]
for i in range(10000):
    result=MarkovChain2 (100,matrix=mat3)
    list2.append(result)
    
print list2
repfreq2=mcStateFreqSum(sims=list2,state=0)
print repfreq2
####The state freq values change, but also seem to converge to a sationary frequency (0.43)

# Now, calculate a vector of probabilities for the focal state (e.g., 'a')
# based on the transition matrix directly (not by simulation). How do these
# values compare to the simulated frequencies?
mat4=[[0.77,0.23],[0.39,0.61]]###This is a variation of my original matrix mat2, just because my Markov Function reverts the propabilities
import numpy
numpy.linalg.matrix_power(M=mat4,n=1000)
#####Also got 0.63 as a probability.

mat5=[[0.23,0.77],[0.61,0.39]]
numpy.linalg.matrix_power(M=mat5,n=1000)
####Got 0.44 as a probability.
####The stationary frequency is quite close to the probabilities calculated.


