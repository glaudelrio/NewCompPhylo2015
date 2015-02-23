# -*- coding: utf-8 -*-"""
"""
Created on Mon Feb 09 23:11:37 2015
@author: Glaucia
"""

import random

mat = [[0.1,0.9],[0.1,0.9]]

import scipy       
def discSamp(events,probs):##importing the discrete sampling function
    ranNum = scipy.random.random()
    cumulProbs = []
    cumulProbs.extend([probs[0]])
    for i in range(1,len(probs)):
        cumulProbs.extend([probs[i]+cumulProbs[-1]])
    for i in range(0,len(probs)):
        if ranNum < cumulProbs[i]:
            return events[i]
    return None
##############################################################################
from numpy.random import uniform
def MarkovChain (x,matrix):##importing the MarkovChain fucntion I made for two states
    rn=uniform(low=0.0, high=1.0, size=None)
    currState=discSamp(events=[0,1],probs=[rn,1-rn])
    list1=[currState]    
    for i in range (x):
        currState=discSamp(events=[currState,1-currState],probs=matrix[currState])
        list1.append(currState)
    return list1

###############################################################################
# ----> Try to finish the above lines before Tues, Feb. 10th <----
# Now try running 100 simulations of 100 steps each. How often does the chain
# end in each state? How does this change as you change the transition matrix?
mat1 = [[0.1,0.9],[0.1,0.9]] ###### 0 and 1 have similar frequencies

sim=[]
for i in range(100):
    r=MarkovChain (100,matrix=mat1)
    p=r[99]
    sim.append(p)
print sim##this is the simulation of MarkovChain 100 times

sim.count(0)##this wiil give the frequencies of ending in each state
sim.count(1)##this wiil give the frequencies of ending in each state

mat2 = [[0.9,0.1],[0.1,0.9]] #### Repetitions of 0 have a large probability here,
sim=[]
for i in range(100):
    r=MarkovChain (100,matrix=mat2)
    p=r[99]
    sim.append(p)
print sim

sim.count(0)##this shows high frequencies of zeros
sim.count(1)

mat3 = [[0.1,0.9],[0.9,0.1]] #### Repetitions of 1 have a large probability here,
sim=[]
for i in range(100):
    r=MarkovChain (100,matrix=mat3)
    p=r[99]
    sim.append(p)
print sim

sim.count(0)
sim.count(1)
###############################################################################
# Try defining a state space for nucleotides: A, C, G, and T. Now define a
# transition matrix with equal probabilities of change between states.
tup=("A","C","T","G")
######Matrix with equal probabilities for each transition
matNucl1=[[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]

######Defining a new function to run Markov Chain using state space with length equal to 4
def MarkovChainNucl (x,matrix,Event):
    tup=tuple(Event)##this will transform a list of events in a tuple
    currState=random.choice(tup)##this will sample a random state from my tuple
    list1=[currState]## creating the list to put the sampled states   
    for i in range (x):##here x is the siz of the chain (you decide x)
        if currState == Event[0]:##this is an if statement to allow choosing the next state according to the respective probabilities in transition matrix
            currState=discSamp(events=Event,probs=matrix[0])##this will sample a new current state according to the respective row in the transition matrix
        elif currState == Event[1]:
            currState=discSamp(events=Event,probs=matrix[1])
        elif currState == Event[2]:
            currState=discSamp(events=Event,probs=matrix[2])
        else:
            currState=discSamp(events=Event,probs=matrix[3])
        list1.append(currState)##Appending all the States in a list
    return list1
#######Testing the function
MarkovChainNucl (10,matrix=matNucl1,Event=["A","C","T","G"])
###############################################################################
# Again, run 100 simulations of 100 steps and look at the ending states. Then
# try changing the transition matrix.
simNucl1=[]
for y in range(100):##simulating for 100 steps
    w=MarkovChainNucl (100,matrix=matNucl1,Event=["A","C","T","G"])
    z=w[99]
    simNucl1.append(z)
print simNucl1

simNucl1.count("A")##Couting the frequencies of the ending nucleotide
simNucl1.count("C")##The frequencies are very similar here, because all the probabilities are the same in the transition matrix
simNucl1.count("T")
simNucl1.count("G")

#######Changing the matrix and getting far more repetitions of C
matNucl2=[[0.25,0.25,0.25,0.25],[0.05,0.85,0.05,0.05],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]##I changed the value in the position Matrix [2][2] to allow repetitions of C

simNucl2=[]
for y in range(100):
    w=MarkovChainNucl (100,matrix=matNucl2,Event=["A","C","T","G"])
    z=w[99]
    simNucl2.append(z)
print simNucl2

#######And also a high frequency of Cs
simNucl2.count("A")
simNucl2.count("C")
simNucl2.count("T")
simNucl2.count("G")

