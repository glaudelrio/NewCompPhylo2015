# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 14:39:59 2015

@author: Glaucia
"""

import scipy as sp  


##Creating a Markov Object with state space a and b 
class MarkovObj(object):
    state=("a","b")##state space
    Probabs=[[0.5,0.5],[0.5,0.5]]##transition matrix (same probabilities for all "changes")
    num=25##size of my Markov Chain

    def dmcSim(self):##this fucntion will simulate a Markov Change with 25 states

    # Define list to hold chain's states
        chain = []  
        def discSamp(events,probs):##discrete sampling function
            ranNum = sp.random.random()
            cumulProbs = []
            cumulProbs.extend([probs[0]])
            for i in range(1,len(probs)):
                cumulProbs.extend([probs[i]+cumulProbs[-1]])
            for i in range(0,len(probs)):
                if ranNum < cumulProbs[i]:
                    return events[i]
        return None
    # Draw a state to initiate the chain
        currState = discSamp(events=self.state,probs=[1.0/len(self.state) for x in self.state])##I used your Markov Chain function because it is my versatile than mine. Mine functions worked only with state spaces (0,1) or (A,C,T,G)
        chain.extend(currState)

    # Simulate the chain over n-1 steps following the initial state
        
        for step in range(1,self.num):
            probs1 = self.Probabs[self.state.index(currState)] # Grabbing row associated with currState
            currState = discSamp(self.state,probs1) # Sample new state
            chain.extend(currState)        
        
        return chain
d=MarkovObj()##assigning a name to the object
print d.state##calling state space inside the object
print d.Probabs##calling the transition matrix
print d.dmcSim()##calling the fucntion with the original attributes.