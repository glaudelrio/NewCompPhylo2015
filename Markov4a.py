# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 19:28:02 2015

@author: Glaucia
"""
"""
In this exercise, we will explore Markov chains that have discrete state spaces
and occur in discrete time steps. To set up a Markov chain, we first need to
define the states that the chain can take over time, known as its state space.
To start, let's restrict ourselves to the case where our chain takes only two
states. We'll call them A and B.
"""
# Create a tuple that contains the names of the chain's states
test=("A","B")
print(test)

"""
The behavior of the chain with respect to these states will be determined by
the probabilities of taking state A or B, given that the chain is currently in
A and B. Remember that these are called conditional probabilities (e.g., the
probability of going to B, given that the chain is currently in state A is
P(B|A).)
We record all of these probabilities in a transition matrix. Each row
of the matrix records the conditional probabilities of moving to the other
states, given that we're in the state associated with that row. In our example
row 1 will be A and row 2 will be B. So, row 1, column 1 is P(A|A); row 1,
column 2 is P(B|A); row 2, column 1 is P(A|B); and row 2, column 2 is P(B|B).
All of the probabilities in a ROW need to sum to 1 (i.e., the total probability
associated with all possibilities for the next step must sum to 1, conditional
on the chain's current state).
In Python, we often store matrices as "lists of lists". So, one list will be
the container for the whole matrix and each element of that list will be
another list corresponding to a row, like this: mat = [[r1c1,r1c2],[r2c1,r2c2]].
We can then access individual elements use two indices in a row. For instance,
mat[0][0] would return r1c1. Using just one index returns the whole row, like
this: mat[0] would return [r1c1,r1c2].
Define a transition matrix for your chain below. For now, keep the probabilties
moderate (between 0.2 and 0.8).
"""
# Define a transition probability matrix for the chain with states A and B
mat = [[0.1,0.9],[0.1,0.9]]
# Try accessing a individual element or an individual row
# Element
####accessing r2c1
mat[1][0]
# Row
####accessing the second row
mat[1]
"""
Now, write a function that simulates the behavior of this chain over n time
steps. To do this, you'll need to return to our earlier exercise on drawing
values from a discrete distribution. You'll need to be able to draw a random
number between 0 and 1 (built in to scipy), then use your discrete sampling
function to draw one of your states based on this random number.
"""
# Import scipy U(0,1) random number generator
from numpy.random import uniform

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
    rn=uniform(low=0.0, high=1.0, size=None)##Getting a random number from an uniform distribution from 0 to 1
    currState=discSamp(events=[0,1],probs=[rn,1-rn])##sampling an initial state between two choices (0,1) and probs are based on the random number sampled on the previous line
    list1=[currState]##creating a list to encompass the states sampled    
    for i in range (x):##this will define the number of steps in your change
        currState=discSamp(events=[currState,1-currState],probs=matrix[currState])##this will sample a new current state according to the respective probabilities in the transition matrix. However, as the order current state changes with the sampling the second row of the transition matrix is reversed.
        list1.append(currState)##appending the states in a list
    return list1##this will return the list of states sampled
##Testing the function with 40 steps
MarkovChain (40,matrix=mat)

# Run a simulation of 10 steps and print the output.
m=MarkovChain (10, matrix=mat)
print m

###############################################################################
