# -*- coding: utf-8 -*-

"""

Created on Sun Mar 08 21:52:37 2015



@author: Glaucia

""" 



from scipy import linalg

import numpy as np



# JMB: Added threshold for diff

def estBrlen (currBr,diff,thresh):

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

    

    

    MargProbsMatrixcurrLike=linalg.expm(A*currBr)##Marginal probability matrix for currBr

    # JMB: I moved this before loop. Matrix exp is only needed once and is expensive.

    for c,m in zip (listfirst,listlast):

        currLikes=MargProbsMatrixcurrLike[c][m]##extracting the corresponding marginal probabilities values of the matrix

        emptylist1.append(currLikes)##putting the probabilities in a list

    currLike=reduce(lambda x, t: x*t,emptylist1)##multipying all the values in the list to get the total Likelihood for the current Branch length

        # JMB: There's no need to multiply the site likelihoods until you've added them all to the list.

        

    while (diff >= thresh):

        MargProbsmatrixupLike=linalg.expm(A*upBr)##doing the same for currBr + diff

        # JMB: Moved before loop

        for c,m in zip (listfirst,listlast):

            upLikes=MargProbsmatrixupLike[c][m]

            emptylist2.append(upLikes)

        upLike=reduce(lambda x, t: x*t,emptylist2)

        # JMB: Same as above - no need to multiply all values until the list is complete.

        MargProbsmatrixdownLike=linalg.expm(A*downBr)##doing the same for currBr - diff   

        for c,m in zip (listfirst,listlast):

            downLikes=MargProbsmatrixdownLike[c][m]

            emptylist3.append(downLikes)

        downLike=reduce(lambda x, t: x*t,emptylist3)

        if (upLike > currLike):

            while (upLike > currLike):

                 currBr=upBr

                 emptylist1=[]

                 emptylist2=[]

                 MargProbsmatrixcurrLike=linalg.expm(A*currBr)

                 for c,m in zip (listfirst,listlast):

                      currLikes=MargProbsmatrixcurrLike[c][m]

                      emptylist1.append(currLikes)

                 currLike=reduce(lambda x, t: x*t,emptylist1)

                 upBr=currBr+diff

                 MargProbsmatrixupLike=linalg.expm(A*upBr)

                 for c,m in zip (listfirst,listlast):

                       upLikes=MargProbsmatrixupLike[c][m]

                       emptylist2.append(upLikes)

                 upLike=reduce(lambda x, t: x*t,emptylist2)

            # return currBr

        elif (downLike > currLike):

            while (downLike > currLike):

                 currBr=downBr 

                 emptylist1=[]

                 emptylist3=[]

                 MargProbsmatrixcurrLike=linalg.expm(A*currBr)

                 for c,m in zip (listfirst,listlast):

                     currLikes=MargProbsmatrixcurrLike[c][m]

                     emptylist1.append(currLikes)

                 currLike=reduce(lambda x, t: x*t,emptylist1)

                 downBr=currBr-diff

                 MargProbsmatrixdownLike=linalg.expm(A*downBr)

                 for c,m in zip (listfirst,listlast):

                     downLikes=MargProbsmatrixdownLike[c][m]

                     emptylist3.append(downLikes)

                 downLike=reduce(lambda x, t: x*t,emptylist3)

            # return currBr

        else:

            "none"

            

        diff *= 0.5 # Reduce diff by 1/2

    return currBr   # Only return value once diff is below thresh

    

import time



estBrlen(currBr=0.3,diff=0.01,thresh=0.01)

estBrlen(currBr=1,diff=0.01,thresh=0.0001)



start=time.clock()

estBrlen(currBr=5,diff=0.01,thresh=0.0000001)

end=time.clock()

print end-start





