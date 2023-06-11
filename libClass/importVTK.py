# ==========================================
# Code created by Leandro Marques at 03/2020
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to import .vtk file


# Converting .msh in a python list

import numpy as np

#Version 2 - Mini element powered
def vtkFile(_file, _numNodes, _numVerts, _numElements, _IEN, _polynomial_option): 
 benchmark_problem = 'Continue simulation'

 #List Python Assembly
 vtkList = [] 
 with open(_file) as vtkFile:
   for line in vtkFile:
    row = line.split()
    vtkList.append(row[:])

 #Mapping 
 #(p.s: How to print index in .txt file python without convert in a list?)
 #To find a string, just make -> in "string" in open(file).read():
 for i in range(0,len(vtkList)):
  if "VECTORS" in vtkList[i]:
    indexVectors = i
    print indexVectors

  if "scalar1" in vtkList[i]:
    indexScalar1 = i
    print indexScalar1

  if "scalar2" in vtkList[i]:
    indexScalar2 = i
    print indexScalar2

  if "scalar3" in vtkList[i]:
    indexScalar3 = i
    print indexScalar3

 
 vx = np.zeros([_numNodes,1], dtype = float)
 vy = np.zeros([_numNodes,1], dtype = float)
 scalar1 = np.zeros([_numNodes,1], dtype = float)
 scalar2 = np.zeros([_numVerts,1], dtype = float)
 scalar3 = np.zeros([_numVerts,1], dtype = float)
 for i in range(0,_numVerts):
  scalar1[i] = float(vtkList[indexScalar1 + 2 + i][0])
  scalar2[i] = float(vtkList[indexScalar2 + 2 + i][0])
  scalar3[i] = float(vtkList[indexScalar3 + 2 + i][0])
  vx[i] = float(vtkList[indexVectors + 1 + i][0])
  vy[i] = float(vtkList[indexVectors + 1 + i][1])


 #Centroid
 if _polynomial_option == 2:
  for e in range(0,_numElements):
   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]
   vx[v4] = (vx[v1] + vx[v2] + vx[v3])/3.0
   vy[v4] = (vy[v1] + vy[v2] + vy[v3])/3.0
   scalar1[v4] = (scalar1[v1] + scalar1[v2] + scalar1[v3])/3.0



 return vx, vy, scalar1, scalar2, scalar3, benchmark_problem


#Version 1 - Loop over list
def vtkFilev1(_file, _polynomial_option): 
 benchmark_problem = 'Continue simulation'

 vtkList = [] 
 with open(_file) as vtkFile:
   for line in vtkFile:
    row = line.split()
    vtkList.append(row[:])

 for i in range(0,len(vtkList)):
  for j in range(0,len(vtkList[i])):
   if vtkList[i][j] == "POINTS":
    numNodes = int(vtkList[i][j+1])

    x = np.zeros([numNodes,1], dtype = float)
    y = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     x[k] = float(vtkList[i+k+1][0])
     y[k] = float(vtkList[i+k+1][1])
    continue  

   if vtkList[i][j] == "CELLS":
    numElements = int(vtkList[i][j+1])

    # Linear Element
    if _polynomial_option == 0 or _polynomial_option == 1:
     IEN = np.zeros([numElements,3], dtype = int)
     polynomial_order = "Linear"
     for e in range(0,numElements):
      IEN[e][0] = int(vtkList[i+e+1][1])
      IEN[e][1] = int(vtkList[i+e+1][2])
      IEN[e][2] = int(vtkList[i+e+1][3])
     continue 

    # Quad Element
    elif _polynomial_option == 3:
     IEN = np.zeros([numElements,6], dtype = int)
     polynomial_order = "Quadratic"
     for e in range(0,numElements):
      IEN[e][0] = int(vtkList[i+e+1][1])
      IEN[e][1] = int(vtkList[i+e+1][2])
      IEN[e][2] = int(vtkList[i+e+1][3])
      IEN[e][3] = int(vtkList[i+e+1][4])
      IEN[e][4] = int(vtkList[i+e+1][5])
      IEN[e][5] = int(vtkList[i+e+1][6])
     continue 


   if vtkList[i][j] == "VECTORS":
    vx = np.zeros([numNodes,1], dtype = float)
    vy = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     vx[k] = float(vtkList[i+k+1][0])
     vy[k] = float(vtkList[i+k+1][1])
    continue  

   if vtkList[i][j] == "scalar1":
    scalar1 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar1[k] = float(vtkList[i+k+2][0])
    continue  

   if vtkList[i][j] == "scalar2":
    scalar2 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar2[k] = float(vtkList[i+k+2][0])
    continue  

   if vtkList[i][j] == "scalar3":
    scalar3 = np.zeros([numNodes,1], dtype = float)
    for k in range(0,numNodes):
     scalar3[k] = float(vtkList[i+k+2][0])
    continue  

 return numNodes, numElements, IEN, x, y, vx, vy, scalar1, scalar2, scalar3, polynomial_order, benchmark_problem
