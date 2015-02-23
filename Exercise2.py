# -*- coding: utf-8 -*-
"""
Created on Sun Jan 25 02:56:50 2015

@author: Glaucia
"""
#*** Discrete Sampling Practice ***
#---> Creating useful functions <---

"""
1.Write a function that multiplies all consecutively decreasing numbers between 
a maximum and a minimum supplied as arguments. (Like a factorial, but not 
necessarily going all the way to 1). This calculation would look like
####
I defined the function MultDecr, which arguments are a maximum and a minimum 
value which express the initial point and the final point of a multiplication 
of a decrescent sequence. After some research, I noticed that a good way to 
write this funtion would be using while loops
"""

def MultDecr(mx,mn):
    n= 1
    while (mx > mn-1 or 0):
        n = n * mx
        mx = mx - 1
    return n

MultDecr(6,2)

"""
2. Using the function you wrote in (1), write a function that calculates the 
binomial coefficient (see Definition 1.4.12 in the probability reading). 
Actually, do this twice. The first time (2a) calculate all factorials fully. 
Now re-write the function and cancel as many terms as possible so you can avoid 
unnecessary multiplication (see the middle expression in Theorem 1.4.13)
####
#2.a
First I defined a function with full factorials 
"""
def Binom(n,k):
   nb1=1
   kb1=1
   sb1=1
   sub = n-k
   while (n > 0):#to calculate the n!
       nb1 = nb1*n
       n = n-1
   while (k > 0):#to calculate the k!
       kb1 = kb1*k
       k= k-1
   while (sub > 0):#to calculate the (n-k)!
       sb1 = sb1*sub
       sub = sub-1 
       result=nb1/(sb1*kb1)# to calculate n!/(n-k)!k!
   return result
        
Binom(5,2)

"""
#2.b
####
And after that I defined a fucntion excluding the multiplications which would 
be canceled by the division
"""

def BinomSimp(n,k):
   nb1=1
   kb1=1
   sub = n-k
   while (n > sub):#here the multiplication goes until n>n-k avoiding some unnecessary multiplications
       nb1 = nb1*n
       n = n-1
   while (k > 0):
       kb1 = kb1*k
       k= k-1
       result=nb1/kb1
   return result
   
BinomSimp(5,2)
"""
#3.Try calculating different binomial coefficients using both the functions from 
(2a) and (2b) for different values of n and k. Try some really big values there 
is a noticeable difference in speed between the (2a) and (2b) function. Which 
one is faster? By roughly how much?
####
As my computer lacks memory, I used not so high values. The fisrt one took 
27 seconds to run and the other took less than two seconds... the second was 
almost 15 times faster than the first one.
"""
teste1=Binom(400,100)
teste2=BinomSimp(400,100)
"""
#4. Use either function (2a) or (2b) to write a function that calculates the 
#probability of k successes in n Bernoulli trials with probability p. This is 
#called the Binomial(n,p) distribution. See Theorem 3.3.5 for the necessary 
#equation. [Hint: pow(x,y) returns x^y (x raised to the power of y).]
####
I used the BinoSimp to calculate the binomial coeficient and then multiplied 
by the variables as presented in theorem 3.3.5
"""
def PMFBinom(n,k,p):
    binomcoef=BinomSimp(n,k)
    return binomcoef*(pow(p,k))*(pow(1-p,n-k))
    
PMFBinom(400,200,0.56)
"""
5. Now write a function to sample from an arbitrary discrete distribution. 
This function should take two arguments. The first is a list of arbitrarily 
labeled events and the second is a list of probabilities associated with these 
events. Obviously, these two lists should be the same length.
####
To write this function I used some previous functions from scipy library. 
the rv_discrete can be used to construct arbitrary discrete distribution based
on the probability of events, since these events are labeled with integer 
numbers.
"""
from scipy.stats import rv_discrete

def DiscreteSample (xk,pk,siz):#I had to use three arguments, since the size of the sample was required
    discrete = rv_discrete(name='discrete', values=(xk, pk))
    sample = discrete.rvs(size=siz)
    x = list(sample)
    return x
####Testing the function
DiceOutcomes=[1,2,3,4,5,6]
Probabilities=[0.167,0.167,0.167,0.167,0.167,0.167] 
DiscreteSample(xk=DiceOutcomes,pk=Probabilities,siz=10)
"""
6. (6) For an alignment of 400 sites, with 200 sites of type 1 and 200 of type
2, sample a new alignment (a new set of site pattern counts) with replacement 
from the original using your function from (5). Print out the counts of the two 
types.
####
Creating a new alignment:
"""
Types=[1,2]#There are to types of sites in the original alignment
Probs=[0.5,0.5]#Each type of site proportion in the alignment represents a probability of 0.5 for each type
newAlignment=DiscreteSample(xk=Types,pk=Probs,siz=400)#with this information I can create a new alignment of 400 sites
print(newAlignment)
newAlignment.count(1)#Counting number of events 1
newAlignment.count(2)#Counting number of events 2
"""
7. Repeat (6) 100 times and store the results in a list:
####
Repeating the sample method 100 times:
"""
list1=[]
for i in range (100):
    outcome=DiscreteSample(xk=[1,2],pk=[0.5,0.5],siz=400)
    counting1=outcome.count(1)
    list1.append(counting1)
print list1
counts100=[float(y) for y in list1]#transforming in floats
####Getting the proportions of number 1 events in 100 alignments ans storing 
####them in a list    
list2=[]
for i in counts100:
    list2.append(i/400)
####Getting the list with the proportions
print list2
"""
8.Of those 100 trials, summarize how often you saw particular proportions of 
type 1 vs. type 2.
#### 
Plotting the frequencies of proportions:
"""
import matplotlib.pyplot as plt
plt.hist(list2)
####Or printing an histogram of the number of successes:
plt.hist(list1)
"""
9.Calculate the probabilities of the proportions you saw in (8) using the 
binomial probability mass function (PMF) from (4).
####
Applying the Probability Mass Function for each of the proportions found on 
the item 9, and storing the results in a list. I considered n=400, k=each 
value in the 100 trials and p equal to 0.5.
"""
list3=[]
for z in list1:
    list3.append(PMFBinom(n=400,k=z,p=0.5))
 
print list3
"""
10. Compare your results from (8) and (9).  
####
As we can see in the histogram, the PMF values (0.01 to 0.01) are much smaller 
than the proportions, also higher values of PMF are more frequent. 
"""
plt.hist(list3)
"""
11.Repeat 7-10, but use 10,000 trials.
####
Doing 10000 trials:
"""
list4=[]
for i in range (10000):
    outcome=DiscreteSample(xk=[1,2],pk=[0.5,0.5],siz=400)
    counting2=outcome.count(1)
    list4.append(counting2)
print list4
plt.hist(list4)#the distribution has "normal shape"
list5=[]
for s in list4:
    list5.append(PMFBinom(n=400,k=s,p=0.5))
print list5
plt.hist(list5)
