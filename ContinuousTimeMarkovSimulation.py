# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 20:11:37 2015
@author: Glaucia
"""
##Continuous time Markov Simulation
import scipy 
import random 
from numpy.random import uniform
from math import log 

class ContMarkSim(object):
    x=10##Ending Time, in this case I chose 10
    Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]]##Q matrix, or Rates matrix
    ##Defining the function that will make my simulation
    def MarkSim(self):
        tup=(0,1,2,3)##Tuple with my state space
        i=random.choice(tup)##Sorting the first state
        E=0##Defining the starting point in the branch
        list1=[]##This will be the list with the waiting times
        list2=[]##This wiil be the list with the probabilities associated with the of-diagonal values
        list3=[]##This will be the list with the possible events for the next change in the chain
        list4=[i]##This will be the list with the states (substitutions)
        def discSamp(events,probs):##defining the function of discrete sample, to sort the events to change
            ranNum = scipy.random.random()
            cumulProbs = []
            cumulProbs.extend([probs[0]])
            for t in range(1,len(probs)):
                cumulProbs.extend([probs[t]+cumulProbs[-1]])
            for t in range(0,len(probs)):
                if ranNum < cumulProbs[t]:
                    return events[t]
            return None
        while E < self.x:##this loop will occur while the sum of the waiting times does not reach the ending time
            u1=uniform(low=0.0, high=1.0, size=None)##sorting an uniform value to give to the next equation
            Wtime=-(1/-(self.Q[i][i]))*log(u1)##calculating the waiting time, by assuming that lambda is equal to the -diagonal (-Qii) value
            list1.append(Wtime)##putting the waiting times in list 1
            E=sum(z for z in list1)##adding the waiting times to make the simulation stop when reaching the ending time
            for g in self.Q[i]:
                if g >= 0:
                    p=-(g/(self.Q[i][i]))##calculating the probabilities of the next events according to the current state. We divide the positive rates of the row, by the diagonal value
                    list2.append(p)##puting the probabilities in a list
            for w in [0,1,2,3]:
                if w != i:##defining the state space of the next change according to the current state
                    list3.append(w)##putting the next state space in a list
            i=discSamp(events=list3,probs=list2)##sampling the next state according to the current state
            list4.append(i)##appending the states in a list
            list5 = [str(h) for h in list4]
            list5 = [z.replace('0', 'A') for z in list5]## Converting the numbers to letters representing nucleotides
            list5 = [a.replace('1', 'C') for a in list5]## Converting the numbers to letters representing nucleotides
            list5 = [u.replace('2', 'G') for u in list5]## Converting the numbers to letters representing nucleotides
            list5 = [n.replace('3', 'T') for n in list5]## Converting the numbers to letters representing nucleotides
        return list5,list1##returning the states(substitutions) and waiting times until reaching an ending time
##Giving a name to the object        
d= ContMarkSim()
##Getting the results of the Continuous time Markov Simultaion
print d.MarkSim()
##Testing changes in the waiting time
d.x=200

