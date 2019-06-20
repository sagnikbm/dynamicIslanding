#10/30/2018 This is the main program for islanding in 200 bus system
# Objective is to maintain load-generation imbalance in each island (and minimize number of lines disconnect)

from gurobipy import *
import numpy as np

##############################################################################################
################################### Load Data ################################################
##############################################################################################
s = open(r'arcs.txt', 'r').read()
arcs, weights= multidict(eval("{" + s + "}"))

s = open(r'genCoherentSteiner.txt', 'r').read()
genCoherentSteiner,genArea= multidict(eval("{" + s + "}"))
#
s = open(r'genSource.txt', 'r').read()
genSource,Area= multidict(eval("{" + s + "}"))

s = open(r'lineSteiner.txt', 'r').read()
lineSteiner,lineArea= multidict(eval("{" + s + "}"))

s = open(r'nodes.txt', 'r').read()
nodes,injection= multidict(eval("{" + s + "}"))

# nodes = range(1,201) # number of nodes
area=range(1,5) # number of areas

# Create optimization model
m = Model('graphPartition')

##############################################################################################
###########################Partition Variables################################################
##############################################################################################

# x(i,h) tells if node i is in area(partition) h
x = m.addVars(nodes, area,vtype=GRB.BINARY, name="x")
# w(i,j,h) is the auxiliary variable indicating if edge (i,j) is in area h
w = m.addVars(arcs,area, vtype=GRB.BINARY, name="w")
# z(i,j) indicates whether edge (i,j) would be added to cost
z = m.addVars(arcs, vtype=GRB.BINARY,name="z")
m.update()

##############################################################################################
###########################Connectivity Variables#############################################
##############################################################################################

# u(j,h) indicates if j_th node is the source node in partition h
u = m.addVars(nodes,area, vtype=GRB.BINARY, name="u")
f1 = m.addVars(arcs,area,name="f")
m.update()


##############################################################################################
##############################Partition Constraints###########################################
##############################################################################################

m.addConstrs(w[i,j,h] <= x[i,h] for i, j in arcs for h in area) #(3a)
m.addConstrs(w[i,j,h] <= x[j,h] for i, j in arcs for h in area) #(3b)
m.addConstrs(z[i,j] <= w.sum(i,j,'*') for i, j in arcs)         #(3d)
m.addConstrs(x.sum(i,'*') == 1 for i in nodes)                  #(3e)
m.addConstrs(x.sum('*',h) >= 6 for h in area)                   #(3f)
m.addConstrs(x[i,j]==1 for i,j in genCoherentSteiner)           #(5)
# m.addConstrs(x.sum('*',h) <= 35 for h in area)                #(3fExtended) ensures maximum of 30 nodes in one area
m.addConstrs(z[i,j]==1 for i,j in lineSteiner)                  #this is an extra constraint added not to cut lines in steiner sub-graph
m.addConstrs(z[i,j] == z[j,i] for i,j in arcs)

# m.addConstrs((quicksum(injection[i]*x[i,h] for i in nodes) <=  80 for h in area)) # gen = load
# m.addConstrs((quicksum(injection[i]*x[i,h] for i in nodes) >= -80 for h in area)) # gen = load

##############################################################################################
##############################Connectivity Constraints########################################
##############################################################################################

m.addConstrs(u[j,h]==1 for j,h in genSource) #(4a-4e) #preselect source node
m.addConstrs(u[j,h]*x.sum('*',h)-x[j,h]+f1.sum('*',j,h)==f1.sum(j,'*',h) for j in nodes for h in area) #connectivity constraint
m.addConstrs(0<=f1[i,j,h] <= 200*z[i,j] for i,j in arcs for h in area) #(4g)
m.addConstrs(f1[i,j,h]>=0 for i,j in arcs for h in area) #(4g)

##############################################################################################
####################################Objective Function########################################
##############################################################################################

#####Add objective as constraints###

t=m.addVars(area, name="t")
m.update()

#main objective function
m.setObjective(t.sum())  #to minimize only load=gen imbalance
#alternative multiple objective function
# m.setObjective(t.sum()+(0.5*(sum(weights[i,j] for i,j in arcs))-z.prod(weights))) #to minimize only load=gen imbalance and no of lines disconnected

#add additional constraints to relax objective function
m.addConstrs((quicksum(injection[i]*x[i,h] for i in nodes) <= t[h] for h in area))
m.addConstrs((quicksum(injection[i]*x[i,h] for i in nodes) >= -t[h] for h in area))


# Write model to text file
m.write('island.lp')


# m.setParam('TimeLimit', 90)

# m.setParam('BestObjStop', 11.17801)

m.modelSense = GRB.MINIMIZE
m.optimize()

#Print solution
print('\nTOTAL COSTS: %g' % m.objVal)
print('SOLUTION:')


##############################Print all variables######################################
for v in m.getVars():
   if v.x==1:
        print('%s %g' % (v.varName, v.x))

##############################Print all variables######################################
   # for v in m.getVars():
   #     print('%s %g' % (v.varName, v.x))

##############################Print particular variables###############################
# for i,j in z:
#         print('Status of Line %s ' % z[i,j])


