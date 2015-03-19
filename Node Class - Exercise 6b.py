"""
Exercise 6 - Creating and Using Node and Tree Classes
@author: jembrown
Below is the beginning of a Node class definition and a simple example of how
to link nodes to form a tree. Use this as a springboard to start thinking about:
- What other attributes of a Node might we like to store?
- How do we define a Tree class? What attributes should it have?
- Can you write a function to print out a parenthetical tree string 
   (e.g., ((spA,spB),spC)) if the only argument passed to the function is a
   root node? This will require recursion.
"""

# ---> Defining Node and Tree classes <---

class Node:
    
    def __init__(self,name="",parent=None,children=None):
        self.name = name
        self.parent = None
        if children is None:
            self.children = []
        else:
            self.children = children
        
        
        
# ---> Creating and linking nodes <---
 
# Creating nodes to build this simple three-taxon tree: ((spA,spB),spC)
       
#  spA     spB  spC
#    \    /     /
#     \  /     /
#      \/     /
#       \    /
#        \  /
#         \/
#         |

# Define the root node to start. It currently has no parents or children.
root = Node("root") 

# Define a node for species C. It is a direct descendant of the root.
spC = Node("SpeciesC",parent=root)
root.children.append(spC)   # Adds spC as a child of the root

# Define a node for the ancestor of species A and B, descending from the root.
ancAB = Node("ancAB",parent=root)
root.children.append(ancAB)
spA = Node("SpeciesA",parent=ancAB) # Creates spA with ancAB as its parent.
spB = Node("SpeciesB",parent=ancAB) # Creates spB with ancAB as its parent.
ancAB.children.append(spA)
ancAB.children.append(spB)


print("ancAB's children: ")
for child in ancAB.children:
    print child.name
    
print("")
print("root's children: ")
for child in root.children:
    print child.name
# Play around with nodes and see if you can build more complicated trees!
###Playing with the Golden-green Woodpecker tree - (P. leucolaemus(P. flavigula(P. aurulentus,P. chrysochloros)))
root = Node("root")
leucolaemus = Node ("leucolaemus",parent=root)
root.children.append(leucolaemus)
ancflav_aur_chry=Node("ancflav_aur_chry",parent=root)
root.children.append(ancflav_aur_chry)
flavigula=Node("flavigula",parent=ancflav_aur_chry)
ancflav_aur_chry.children.append(flavigula)
ancaur_chry=Node("ancaur_chry",parent=ancflav_aur_chry)
ancflav_aur_chry.children.append(ancaur_chry)
aurulentus=Node("aurulentus",parent=ancaur_chry)
ancaur_chry.children.append(aurulentus)
chrysochloros=Node("chrysochloros",parent=ancaur_chry)
ancaur_chry.children.append(chrysochloros)

print("ancaur_chry children: ")
for child in ancaur_chry.children:
    print child.name
print("root's children: ")
for child in root.children:
    print child.name
print("ancflav_aur_chry children: ")
for child in ancflav_aur_chry.children:
    print child.name  

# Eventually, we will want to create a Tree class, where a parenthetical tree
# string is passed as an argument to the constructor and it automatically creates
# all the nodes and links them together. Start thinking about how to do that.


# Let's go ahead and define a Tree object that houses all these nodes and 
# organizes methods associated with them.
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
    
    def __init__ (self,NodeChain=1,v=10,Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]],waittimes=[],list4=[], states=[],TotalProbs=[],MargProb=[],Likelihood=[]):
        self.NodeChain=NodeChain
        self.v=v##Ending Time, in this case I chose 10
        self.Q=Q##Q matrix, or Rates matrix
        ##Tuple with my state space
        self.waittimes=[]##This will be the list with the waiting times
        self.list4=[]##This will be the list with the states (in numbers where 0=A,1=C,2=G and 3=T)
        self.states=states##This will give you the states as nucleotides.
        self.TotalProbs=TotalProbs
        self.MargProb=MargProb
        self.Likelihood=Likelihood
    ##Defining the function that will make my simulation
    def MarkSim(self):
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
        list6=[]
        list7=[]
        list10=[]
        E=0##Defining the starting point in the branch        
        i=self.NodeChain
        self.list4=[i]
        M4=StationaryProbs[1][i]
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




# We'll need this later, I promise. It's always better to put import statements
# outside of class definitions. In fact, it's best to put them all at the top
# of a file. This imports the ctmc class that we previously defined.

class Tree:
    """
    Defines a class of phylogenetic tree, consisting of linked Node objects.
    """
    
    def __init__(self):
        """
        The constructor really needs to be more flexible, but for now we're 
        going to define the whole tree structure by hand. This just uses
        the same statements we used above. By next Thurs (3/19), see if you can
        write a constructor that takes a parenthetical tree as its argument and 
        builds the corresponding tree in memory. 
        """
        self.root = Node("root") 
        self.spC = Node("SpeciesC",parent=self.root)
        self.root.children.append(self.spC)
        self.ancAB = Node("ancAB",parent=self.root)
        self.root.children.append(self.ancAB)
        self.spA = Node("SpeciesA",parent=self.ancAB)
        self.spB = Node("SpeciesB",parent=self.ancAB)
        self.ancAB.children.append(self.spA)
        self.ancAB.children.append(self.spB)
        # Now, let's add branch lengths to our Node objects (remember, these fields
        # can be added arbitrarily in Python). In the future, we should probably include
        # branch lengths in the Node constructor.
        self.spA.brl = 0.1
        self.spB.brl = 0.1
        self.spC.brl = 0.2
        self.ancAB.brl = 0.1
        self.root.brl = 0
        # We're also going to add lists to each node that will hold simulated
        # sequences.
        self.spA.seq = []
        self.spB.seq = []
        self.spC.seq = []
        self.ancAB.seq = []
        self.root.seq = []


    # Write a recursive function that takes the root node as its only argument and
    # prints out all the names of the terminal nodes in the tree. Due next Tues (3/17).
    def printNames(self,node):
        """
        A method of a Tree object that will print out the names of its
        terminal nodes.
        """
        if len(node.children) > 0:
            for child in node.children:
                self.printNames(child)
        else:
            print node.name
    # Write a recursive function to calculate the total tree length (the sum of
    # all the branch lengths). Again, the root node of a tree should be the only 
    # argument the first time this function is called. Due next Tues (3/17).
    def treeLength(self,node):
        """
        A method to calculate and return total tree length.
        """
        sumbrl=0
        if node.children is not None:
            for child in node.children:
                sumbrl += self.treeLength(child)###this sums the value of branch length and assigns the new value to the sum
        return node.brl + sumbrl
 

    # Write a recursive function that takes the root node as one of its arguments
    # and prints out a parenthetical (Newick) tree string. Due next Tues (3/17).
    
    def newickLists(self,node):##this function gives a tree with relationships in a list of lists
        """
        A method of a Tree object that will print out the Tree as a 
        parenthetical List.
        """
        
        listoflists=[]
        if len(node.children) > 0:
            for child in node.children:
                d=self.newickLists(child)
                listoflists.append(d)
        else:
            listoflists.append(node.name)
        return listoflists###instead of a string this returns a list of lists, I know this could be a problem in the future, but for now the relationships of the nodes are mantained
   
   ##This is the recursive function that will print the parenthetical tree in a string
   ##Based on Subir's code
    
    def newick(self,node):
        """
        A method of a Tree object that will print out the Tree as a 
        parenthetical string (Newick format).
        """
        parenthetical = "(" 
        if len(node.children) == 0:
            return node.name + ":" + str(node.brl)
        else:
            for child in node.children:
                if node.children[-1] == child: 
                    parenthetical += self.newick(child)
                else:
                    parenthetical += self.newick(child) + ","
            if node.brl is not 0:
                parenthetical += "):" + str(node.brl)
            else:
                parenthetical += ")"
            return parenthetical    



    # Now, let's write a recursive function to simulate sequence evolution along a
    # tree. This amounts to simply simulating evolution along each branch 
    # from the root towards the tips. We'll need to use our ctmc class for setting the 
    # conditions of our simulation, which is why we imported it above our tree 
    # class definition. In this case, we've stored the definition of our ctmc 
    # class in a separate file (ctmc.py) to keep our tree code compact.
    # Now, let's add a ctmc object to each internal node in our tree (except the
    # root). Again, it would be best to add the ctmcs as part of the Node
    # constructor, if we know that we'll be simulating data.
    
    # Try to get this simulator and associated functions working by next Thurs. (3/19)    
        """
        #I did not use this function. My simulation in a tree is in the next two functions
        This method of a Tree object defines a ctmc object associated with all
        nodes that have a branch length (i.e., all but the root).

        def setModels(self,node):
      
        if node.children is not None:
            for child in node.children:
                if child.brl > 0:
                    ContMark=CTMC()
                return ContMark
        """
    def RootState(self):##This function draws a starting state from the stationary frequencies.
        Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]]
        A = np.squeeze(np.asarray(Q))
        StationaryProbs=linalg.expm(A*1000)##getting the matrix with stationary probability. 
        StartingPointvalue=discSamp(events=[0,1,2,3],probs=StationaryProbs[1])
        return StartingPointvalue## 0 corresponds to A, 1 corresponds to C, 2 corresponds to G and 3 to T
    
    def simulate(self,node,StartingPoint):##This function simulates markov chains in each node of my tree
        """
        This method simulates evolution along the branches of a tree, taking
        the root node as its initial argument.
        """
        if node.parent is None:##dealing with the root

            for child in node.children: 
                chain=ContMarkSim(v=child.brl,NodeChain=StartingPoint)##the first starting point will come from the RootState function        
                chain.MarkSim()
                node.chain = chain.list4[-1]##this gives me the last nucleotide  of the chain
                print chain.list4
                self.simulate(node=child,StartingPoint=node.chain)##recursion starting from the last state of the previous node
        else:
            #if len(node.children)>0:
            for child in node.children:##dealing with internal nodes
                chain=ContMarkSim(v=child.brl,NodeChain=node.parent.chain)##the first state will be the last one of the parent chain    
                chain.MarkSim()
                node.chain = chain.list4[-1]
                print chain.list4
                self.simulate(node=child,StartingPoint=node.chain)
             
    
    #def printSeqs(self,node): My last function prints the sequences for each node.
        """
        This method prints out the names of the tips and their associated
        sequences as an alignment (matrix).
        """
d=Tree()
d.printNames(node=d.root)
d.treeLength(node=d.root)
d.newickLists(node=d.root)##This gives a list of lists and not a string...=P
d.newick(node=d.root)
d.RootState()
d.simulate(node=d.root,StartingPoint=d.RootState())##printing the sequences for each node.