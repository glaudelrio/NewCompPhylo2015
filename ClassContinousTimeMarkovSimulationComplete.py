# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 08:05:02 2015

@author: Glaucia
"""

##Continuous time Markov Simulation
import scipy  
from numpy.random import uniform
from math import log 
from math import exp
from itertools import tee, islice, chain, izip
import operator
import functools
from scipy import linalg
import numpy as np



class ContMarkSim(object):
    
    def __init__ (self,v=10,Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]],waittimes=[],list4=[], states=[],TotalProbs=[],MargProb=[]):
        self.v=v##Ending Time, in this case I chose 10
        self.Q=Q##Q matrix, or Rates matrix
        ##Tuple with my state space
        self.waittimes=waittimes##This will be the list with the waiting times
        self.list4=list4##This will be the list with the states (in numbers where 0=A,1=C,2=G and 3=T)
        self.states=states##This will give you the states as nucleotides.
        self.TotalProbs=TotalProbs
        self.MargProb=MargProb
    ##Defining the function that will make my simulation
    def MarkSim(self):
        statespace=[0,1,2,3]
        ##defining the function of discrete sample, to sort the events to change
        def discSamp(events,probs):
            ranNum = scipy.random.random()
            cumulProbs = []
            cumulProbs.extend([probs[0]])
            for t in range(1,len(probs)):
                cumulProbs.extend([probs[t]+cumulProbs[-1]])
            for t in range(0,len(probs)):
                if ranNum < cumulProbs[t]:
                    return events[t]
            return None
        A = np.squeeze(np.asarray(self.Q))##this transforms my Q matrix in an array A
        StationaryProbs=linalg.expm(A*1000)##getting the matrix with stationary probability. 
        i=discSamp(events=statespace,probs=StationaryProbs[1])##Sorting the first state
        E=0##Defining the starting point in the branch
        self.list4=[i]
        list6=[]
        list7=[]
        list10=[] 
        M4=StationaryProbs[1][i]##This will give the respective stationary probability to sample the first state of my chain
        ##this for loop will calculate the transition matrix for each state from the Q matrix
        for d in range(len(self.Q)): 
            list11=[]    
            for g in self.Q[d]:
                if g > 0: 
                    p=-(g/(self.Q[d][d]))
                else:
                    p = 0
                list11.append(p)
        list10.append(list11)
        T= list10##this is the transition matrix
        ##this fucntion will be used to set two lists, one with the current sampled states, and another with the respective next states. This will be used to calculate the probabilities of each state change.
        def previous_and_next(some_iterable):
            prevs, items, nexts = tee(some_iterable, 3)
            prevs = chain([None], prevs)
            nexts = chain(islice(nexts, 1, 0), [0])
            return izip(prevs, items, nexts)
        while E < self.v:##this loop will occur while the sum of the waiting times does not reach the ending time
            list2=[] 
            u1=uniform(low=0.0, high=1.0, size=None)##sorting an uniform value to give to the next equation
            Wtime=-(1/-(self.Q[i][i]))*log(u1)##calculating the waiting time, by assuming that lambda is equal to the -diagonal (-Qii) value
            self.waittimes.append(Wtime)##putting the waiting times in list 1
            E=sum(z for z in self.waittimes)##adding the waiting times to make the simulation stop when reaching the ending time
            for g in self.Q[i]:
                if g > 0: 
                    p=-(g/(self.Q[i][i]))##calculating the probabilities of the next events according to the current state. We divide the positive rates of the row, by the diagonal value
                else:
                    p = 0
                list2.append(p)##puting the probabilities in a list
            for s in self.waittimes:
                exponPDF = -(self.Q[i][i]) * exp(- (-(self.Q[i][i]))*s)##this function calculates the exponential probability density function for each waiting time according to the current state
                list6.append(exponPDF)##this list keeps exponPDF values
                M1=functools.reduce(operator.mul, list6, 1)## this list has the multiplication of all exponPDF values
            i=discSamp(events=[0,1,2,3],probs=list2)##sampling the next state according to the current state
            self.list4.append(i)##appending the states in a list
            list67=[]
            list68=[]
            for previous, item, nxt in previous_and_next(self.list4):
                list67.append(item)##this list will keep the previous sampled nucleotide
                list68.append(nxt)##this list will keep the subsequent nucleotide
            for c,m in zip (list67,list68):##this loop will give the right values to transition matrix values
                if m != list68[-1]:
                    statesprobs=T[c][m]##this will give the probabilities of going from the previous state to the next
                    list7.append(statesprobs)##this will append the probabilities of each state in list7
                else:
                    pass
                    M2=functools.reduce(operator.mul, list7, 1)#this will multiply the probabilities of sampling each state
            MargProbs=linalg.expm(A*self.v)##this calculates the marginal probabilities
            finaltime=self.waittimes[-1]
            finalstate=self.list4[-1]
            cdffinaltime=1-(2.71826**(-(self.Q[finalstate][finalstate])*finaltime))##calculating the probability of last waiting time without a particular change in it.
            M3=1-cdffinaltime
            self.MargProb=MargProbs[self.list4[0]][self.list4[-1]]
            self.TotalProbs=M4*M1*M2*M3##this multiplies the probabilities for each waiting time and for the each nucleotide, in other words, the probability of all the events that occurred in my simulation
            self.states = [str(h) for h in self.list4]
            self.states = [z.replace('0', 'A') for z in self.states]## Converting the numbers to letters representing nucleotides
            self.states = [a.replace('1', 'C') for a in self.states]## Converting the numbers to letters representing nucleotides
            self.states = [u.replace('2', 'G') for u in self.states]## Converting the numbers to letters representing nucleotides
            self.states = [n.replace('3', 'T') for n in self.states]## Converting the numbers to letters representing nucleotides
        return self.states,self.waittimes,self.TotalProbs,self.MargProb##returning the states(substitutions) and waiting times until reaching an ending time




##Giving a name to the object        
d= ContMarkSim()
##Getting the results of the Continuous time Markov Simultaion
print d.MarkSim()
print d.states
print d.waittimes
print d.TotalProbs
print d.MargProb
##Testing changes in the waiting time
banana=ContMarkSim(v=200)
print banana.MarkSim()
print banana.states
print d.waittimes
print d.TotalProbs
print d.MargProb


####trying to test a series of possible branch lengths to go from a nucleotide to another... I put this outside the method in order to think in a way of working with a character...        
Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]] 
A = np.squeeze(np.asarray(Q))##transforming the Q matrix in an array 
Liks=[]##list for likelihoods  
for y in range(100):##varying the branch length from 0 to 99
    MargProbsmatrix=linalg.expm(A*y)##getting the matrix with marginal probabilities
    probsmargAG_CT_AT_AA_TT_CC_GG=(MargProbsmatrix[0][2])*(MargProbsmatrix[1][3])*(MargProbsmatrix[0][3])*(MargProbsmatrix[0][0])*(MargProbsmatrix[3][3])*(MargProbsmatrix[1][1])*(MargProbsmatrix[2][2])##calculating the probabilities of starting at one nucleotide and going to other for 4 sites
    Liks.append(probsmargAG_CT_AT_AA_TT_CC_GG)  
print Liks
Liks.index(max(Liks))##the best branch length (with maximum likelihood) is the biggest one (v=99) even when using 4 nucleotides change


###########I can make a function to get the maximum likelihood estimator
def mle(vCurr,diff,Q,states):##
    first=states[0]
    last=states[-1]    
    A = np.squeeze(np.asarray(Q))    
    vUp=vCurr+diff
    vDown=vCurr-diff
    MargProbsmatrixvCurr=linalg.expm(A*vCurr)
    LvCurr=MargProbsmatrixvCurr[first][last]
    MargProbsmatrixvUp=linalg.expm(A*vUp)
    LvUp=MargProbsmatrixvUp[first][last]
    MargProbsmatrixvDown=linalg.expm(A*vDown)
    LvDown=MargProbsmatrixvDown[first][last]
    if (LvCurr < LvUp):    
        while (LvCurr < LvUp):
            vCurr=vUp
            MargProbsmatrixvCurr=linalg.expm(A*vCurr)
            LvCurr=MargProbsmatrixvCurr[first][last]
            vUp=vCurr+diff
            MargProbsmatrixvUp=linalg.expm(A*vUp)
            LvUp=MargProbsmatrixvUp[first][last]
        return vCurr
    else:
        while (LvCurr < LvDown):
            vCurr=vDown
            MargProbsmatrixvCurr=linalg.expm(A*vCurr)
            LvCurr=MargProbsmatrixvCurr[first][last]
            vDown=vCurr-diff
            MargProbsmatrixvDown=linalg.expm(A*vDown)
            LvDown=MargProbsmatrixvDown[first][last]
        return vCurr
mle(4,1,[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]],[0,2,1,3,1,0,1,2] )