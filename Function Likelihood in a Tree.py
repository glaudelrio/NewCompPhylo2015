# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 19:03:59 2015

@author: Glaucia
"""

from scipy import linalg
import numpy as np
import re

class Node:
    
    def __init__(self,name="",parent=None,children=None):
        self.name = name
        self.parent = None
        if children is None:
            self.children = []
        else:
            self.children = children

extVar = None
###This function reads a newick and determines the relationship among nodes. If you want more details or comments, go to the Node Class - Exercise6b    
def ReadNewick (data,base,space = None, extVar = None):

        if base == "root" :#the content of this if statement will deal with the sons of the root
            root = Node("root")
            root.brl = "0"
            extVar = root
            son1 = data.partition('(')[-1].rpartition(')')[0]
            regEx1 = re.compile(r'(.*?)\(.*\)')
            result1 = re.findall(regEx1, son1)
            result1=result1.pop(0)
            result2 = son1.replace(result1,"")
            result1 = result1[:-1]
            result1=result1.split(",") 
            result1.append(result2)
            for obj in result1:
                branchname=obj.rpartition(":")[0]
                child=Node(name=branchname,parent="root")
                print child.name #printing the sons of the root
                root.children.append(child)
                brl=obj.rpartition(":")[-1]
                child.brl=brl
                print child.brl #and their brancj length
                child.parent = root
                print child.parent #parent
            for item in root.children:
                if item.name[0]=="(":
                    ReadNewick(data=item.name,base="x",space=item)
            return extVar
        else:
            son2 = data.partition('(')[-1].rpartition(')')[0]
            regEx1 = re.compile(r'(.*?)\(.*\)')
            result3 = re.findall(regEx1, son2)
            if result3 != [] and result3 != ['']:
                result3=result3.pop(0)
                result4 = son2.replace(result3,"")
                result3=result3.split(",")
                result3 = result3[:-1]
                result3.append(result4)
                for obj in result3:
                    branchname=obj.rpartition(":")[0]
                    child=Node(name=branchname,parent=data)
                    print child.name #this will give the child of each internal node
                    space.children.append(child)
                    brl=obj.rpartition(":")[-1]
                    child.brl=brl
                    print child.brl #this will give the branch length for each internal node
                    child.parent=space
                    print child.parent #parent
                    if child.name[0]=="(":
                        ReadNewick(data=child.name,base="x",space=child)
            else:
                if son2.count("(") > 0:
                    son3=son2.split(",(")
                    print son3
                    for value in son3:
                        if value[0]!="(":
                            value1 = son3.pop(son3.index(value))
                            value2 = "(" + value1
                            son3.append(value2)
                    print son3
                    for obj in son3:
                        branchname=obj.rpartition(":")[0]
                        child=Node(name=branchname,parent=data)
                        space.children.append(child)
                        brl=obj.rpartition(":")[-1]
                        child.brl=brl
                        print child.name
                        print child.brl
                        child.parent=space 
                        print child.parent #parent
                        ReadNewick(data=child.name,base="x",space=child)
                else:
                    son3 = son2.split(",")
                    print son3,"teste"
                    for obj in son3:
                        branchname=obj.rpartition(":")[0]
                        child=Node(name=branchname,parent=data)
                        space.children.append(child)
                        brl=obj.rpartition(":")[-1]
                        child.brl=brl
                        print child.name
                        print child.brl
                        child.parent=space
                        print child.parent #parent
                    
                    
SomeTree="(A:1,(B:0.1,(C:0.2,D:1):1):0.6)" 
x=ReadNewick(data=SomeTree,base="root",space=None)

##just a function to associate some nucleotides with the tips of a tree, this function is not flexible. I coul make it more flexible and able to deal with other kinds of trees and sequences, but I will leave it like this for now just to help to calculate the likelihood 
def printNucl(node):
        """
        A method of a Tree object that will print out the names of its
        terminal nodes.
        """
        if len(node.children) > 0:
            node.nucl = []
            for child in node.children:
                printNucl(child)
        else:
            if node.name == "A":
                node.nucl=[1,0,0,0]#terminal A has nucleotide A
            elif node.name == "B":
                node.nucl=[1,0,0,0]#terminal B has nucleotide A
            elif node.name == "C":
                node.nucl=[1,0,0,0]#terminal C has nucleotide A
            elif node.name == "D":
                node.nucl=[1,0,0,0]#terminal D has nucleotide A
            print node.nucl
        return node
                
            

x = printNucl(x)
######################################
#####This function calculates de likelihood of a Tree for one nucleotide, starting with just the root and a Q matrix
#####If you want to deal with the likelihood for several nucleotides, you can multipliy the results of likelihood obtained in each tree 
######################################
def TreeLikelihood (node,Q):
    array = np.squeeze(np.asarray(Q))##this transforms my Q matrix in an array A    
    if len(node.children) > 0 and node.nucl == []:##I made this if statement to get from the root to the tips of the tree. node.nucl represents a list of 4 values, where each value is the probabilitie of having a nucleotide at that site.
       for child in node.children:
           TreeLikelihood(child,Q)#listA,listC,listG,listT)
    
    else:
        listA=[]#this list will keep the information on the probabilities of having an A at each internal node
        listC=[]#this list will keep the information on the probabilities of having an C at each internal node
        listG=[]#this list will keep the information on the probabilities of having an G at each internal node
        listT=[]#this list will keep the information on the probabilities of having an T at each internal node
        if node.name != "root":#to deal with all the other nodes except for the root
            for x in node.parent.children:#dealing with the children of the node's parent, in other words, the sisters of the node
                print "node name:", x.name#just checking the node
                print "node branch length:",x.brl,#checking its branch length
                print "probabilities of each nucleotide for this node:",x.nucl#checking the presence of that list with the probabilities for each nucleotide
                if x.nucl != []:#if the node has some information on the probabilities of each nucleotide do the folowing calculations:
                    MP=linalg.expm(array*float(x.brl))#getting the marginal probabilities in accordance with the node branch length
                    x.resA=MP[0][0]*x.nucl[0]+MP[0][1]*x.nucl[1]+MP[0][2]*x.nucl[2]+MP[0][3]*x.nucl[3]#Paa*Pa+Pac*Pc+Pag*Pg+Pat*Pt
                    x.resC=MP[1][0]*x.nucl[0]+MP[1][1]*x.nucl[1]+MP[1][2]*x.nucl[2]+MP[1][3]*x.nucl[3]#Pca*Pa+Pcc*Pc+Pcg*Pg+Pct*Pt
                    x.resG=MP[2][0]*x.nucl[0]+MP[2][1]*x.nucl[1]+MP[2][2]*x.nucl[2]+MP[2][3]*x.nucl[3]#Pga*Pa+Pgc*Pc+Pgg*Pg+Pgt*Pt
                    x.resT=MP[3][0]*x.nucl[0]+MP[3][1]*x.nucl[1]+MP[3][2]*x.nucl[2]+MP[3][3]*x.nucl[3]#Pta*Pa+Ptc*Pc+Ptg*Pg+Ptt*Pt
                    listA.append(x.resA)#putting the results in a list for A
                    listC.append(x.resC)#putting the results in a list for C
                    listG.append(x.resG)#puttinh the results in a list for G
                    listT.append(x.resT)#putting the results in a list for T
                if len(listA)>1:##After doing for both sisters, the list will have more than one result.
                    squareA=reduce(lambda z, y: z*y,listA)#so we just need to multipliy the values inside the list to have the probabilities for the internal node.
                if len(listC)>1:##doing this for all the nucleotides
                    squareC=reduce(lambda z, y: z*y,listC)
                if len(listG)>1:                
                    squareG=reduce(lambda z, y: z*y,listG)
                if len(listT)>1:
                    squareT=reduce(lambda z, y: z*y,listT)
                    lista = [squareA,squareC,squareG, squareT]
                    x.parent.nucl=lista##than we can associate the probabilities for each nucleotide to the parent of the node (an internal node)
                    print "probability for each nucleotide associated with the node's parent:",x.parent.nucl
                    if x.parent.nucl != []:#now that the parent has some probabilities associated with it we just need to do the recursion now using the parent node
                            TreeLikelihood (node.parent,Q)#recursion with the parent node
                    else:
                        pass 
      
        else:#to deal with the root
            StatProb=linalg.expm(array*100)#we just need to calculate the stationary probabilities using a high value as branch length
            likelihood=node.nucl[0]*StatProb[0][0]+node.nucl[1]*StatProb[0][1]+node.nucl[2]*StatProb[0][2]+node.nucl[3]*StatProb[0][3]#and multipliy the probabilities associated with the root by the respective stationary probability and sum all the values obtained
            print "LIKELIHOOD:",likelihood #this will give the likelihood          
##As all the terminals have A as a nucleotide the last list (before multiplying the values by the stationary probs), you can see that the value corresponding to the probability of having an A in the root is much bigger than the other values: [0.0031342687620121074, 0.00031325924070299733, 0.00047299874512964405, 0.00093586514844217763]        
TreeLikelihood(node=x,Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]])
