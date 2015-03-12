# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 20:39:05 2015

@author: Glaucia
"""

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
        
banana=Node()
banana.name="picles"      
banana.name        
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
spC = Node("Species C",parent=root)
root.children.append(spC)   # Adds spC as a child of the root

# Define a node for the ancestor of species A and B, descending from the root.
ancAB = Node("ancAB",parent=root)
root.children.append(ancAB)
spA = Node("Species A",parent=ancAB) # Creates spA with ancAB as its parent.
spB = Node("Species B",parent=ancAB) # Creates spB with ancAB as its parent.
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
#(leucolaemus(flavigula(aurulentus,chrysochloros)))
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
##Examples from the book:
tax_dict = {
'Pan troglodytes' : 'Hominoidea',
'Pongo abelii' : 'Hominoidea',
'Hominoidea' : 'Simiiformes',
'Simiiformes' : 'Haplorrhini',
'Tarsius tarsier' : 'Tarsiiformes',
'Haplorrhini' : 'Primates',
'Tarsiiformes' : 'Haplorrhini',
'Loris tardigradus' : 'Lorisidae',
'Lorisidae' : 'Strepsirrhini',
'Strepsirrhini' : 'Primates',
'Allocebus trichotis' : 'Lemuriformes', 'Lemuriformes' : 'Strepsirrhini', 'Galago alleni' : 'Lorisiformes',
'Lorisiformes' : 'Strepsirrhini',
'Galago moholi' : ' Lorisiformes' }

def get_ancestors(taxon):
    first_parent = tax_dict.get(taxon) 
    second_parent = tax_dict.get(first_parent) 
    third_parent = tax_dict.get(second_parent) 
    return[first_parent, second_parent, third_parent]
    
def get_ancestors2(taxon): 
    result = [taxon]
    while taxon != 'Primates':
        result.append(tax_dict.get(taxon)) 
        taxon = tax_dict.get(taxon) 
    return result

def get_ancestors3(taxon): 
    if taxon == 'Primates': 
        return [taxon]
    else:
        parent = tax_dict.get(taxon) 
        parent_ancestors = get_ancestors(parent) 
        return [parent] + parent_ancestors
        
def get_ancestors4(taxon):
    print('calculating ancestors for ' + taxon) 
    if taxon == 'Primates':
        print('taxon is Primates, returning an empty list') 
        return []
    else:
        print('taxon is not Primates, looking up the parent') 
        parent = tax_dict.get(taxon) 
        print('the parent is ' + parent + ' ') 
        print('looking up ancestors for ' + parent) 
        parent_ancestors = get_ancestors(parent) 
        print('parent ancestors are ' + str(parent_ancestors)) 
        result = [parent] + parent_ancestors 
        print('about to return the result: ' + str(result)) 
        return result

get_ancestors('Galago alleni')

######################################################################################################################
##Trying to find a way to transform parenthetical trees in a dictionary. I found some piece of code in Stack Overflow
def toTree(expression):
    tree = dict()
    msg =""
    stack = list()
    for char in expression:
        if(char == '('):
            stack.append(msg)
            msg = ""
        elif char == ')':
            parent = stack.pop()
            if parent not in tree:
                tree[parent] = list()
            tree[parent].append(msg)
            msg = parent
        else:
            msg += char
    return tree

expression = "(Root(AB(ABC)(CBA))(CD(CDE)(FGH)))"
print toTree(expression)

##Now I just need to join this function with one of the fucntions of the Advanced Python for Biologists Chapter 2 in a Method.

