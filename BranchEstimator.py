# -*- coding: utf-8 -*-
"""
Created on Sun Mar 08 21:52:37 2015

@author: Glaucia
""" 

from scipy import linalg
import numpy as np

def estBrlen (currBr,diff):
    upBr=currBr+diff
    downBr=currBr-diff    
    sites=[[0,1,3,1,2,0],[1,2,1,3,2,3,0,2,0],[0,3,2,3,1],[1,0,2,1],[3,2,1,3],[2,1,2,3,1,2],[1,3,2,0,2,1],[2,3,1,0,2,0],[0,1,2,3,0,1]]##chains of nucleotide states
    Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]] ##Q matrix
    A = np.squeeze(np.asarray(Q))##Q matrix as array
    listfirst=[]
    listlast=[]    
    for y in sites:
        first=y[0]##getting the fisrt value of each site
        last=y[-1]
        listfirst.append(first)
        listlast.append(last)
    emptylist1=[]
    emptylist2=[]
    emptylist3=[]
    for c,m in zip (listfirst,listlast):
        MargProbsMatrixcurrLike=linalg.expm(A*currBr)##Marginal probability matrix for currBr
        currLikes=MargProbsMatrixcurrLike[c][m]##extracting the corresponding marginal probabilities values of the matrix
        emptylist1.append(currLikes)##putting the probabilities in a list
        currLike=reduce(lambda x, t: x*t,emptylist1)##multipying all the values in the list to get the total Likelihood for the current Branch length
    for c,m in zip (listfirst,listlast):
        MargProbsmatrixupLike=linalg.expm(A*upBr)##doing the same for currBr + diff
        upLikes=MargProbsmatrixupLike[c][m]
        emptylist2.append(upLikes)
        upLike=reduce(lambda x, t: x*t,emptylist2)
    for c,m in zip (listfirst,listlast):
        MargProbsmatrixdownLike=linalg.expm(A*downBr)##doing the same for currBr - diff
        downLikes=MargProbsmatrixdownLike[c][m]
        emptylist3.append(downLikes)
        downLike=reduce(lambda x, t: x*t,emptylist3)      
    if (upLike > currLike):
        while (upLike > currLike):
             currBr=upBr
             emptylist1=[]
             emptylist2=[]
             for c,m in zip (listfirst,listlast):
                  MargProbsmatrixcurrLike=linalg.expm(A*currBr)
                  currLikes=MargProbsmatrixcurrLike[c][m]
                  emptylist1.append(currLikes)
                  currLike=reduce(lambda x, t: x*t,emptylist1)
             for c,m in zip (listfirst,listlast):
                   upBr=currBr+diff
                   MargProbsmatrixupLike=linalg.expm(A*upBr)
                   upLikes=MargProbsmatrixupLike[c][m]
                   emptylist2.append(upLikes)
                   upLike=reduce(lambda x, t: x*t,emptylist2)
        return currBr
    elif (downLike > currLike):
        while (downLike > currLike):
             currBr=downBr 
             emptylist1=[]
             emptylist3=[]
             for c,m in zip (listfirst,listlast):
                 MargProbsmatrixcurrLike=linalg.expm(A*currBr)
                 currLikes=MargProbsmatrixcurrLike[c][m]
                 emptylist1.append(currLikes)
                 currLike=reduce(lambda x, t: x*t,emptylist1)
             for c,m in zip (listfirst,listlast):
                 downBr=currBr-diff
                 MargProbsmatrixdownLike=linalg.expm(A*downBr)
                 downLikes=MargProbsmatrixdownLike[c][m]
                 emptylist3.append(downLikes)
                 downLike=reduce(lambda x, t: x*t,emptylist3)
        return currBr
    else:
        "none"
    
    
estBrlen(currBr=0.3,diff=0.01)
estBrlen(currBr=1,diff=0.01)
estBrlen(currBr=5,diff=0.01)