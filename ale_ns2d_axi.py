# =======================
# Importing the libraries
# =======================

import os
initial_path = os.getcwd()

import sys
folderClass = './libClass'
sys.path.insert(0, folderClass)

from tqdm import tqdm
from time import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

import searchMSH
import importMSH
import assembly
import benchmarkProblems
import importVTK
import ALE	
import semiLagrangian
import exportVTK
import relatory



print '''
                 COPYRIGHT                    
 ===========================================
 Simulator: %s
 created by Leandro Marques at 06/2023
 e-mail: leandro.marques@mecanica.coppe.ufrj.br
 COPPE/Departament of Mechanical Engineering
 Federal University of the Rio de Janeiro
 ===========================================
\n''' %sys.argv[0]



print ' ------'
print ' INPUT:'
print ' ------'

print ""


print ' ----------------------------------------------------------------------------'
print ' (0) - Debug'
print ' (1) - Simulation'
simulation_option = int(raw_input("\n enter simulation option above: "))
print' ----------------------------------------------------------------------------\n'



print ' ----------------------------------------------------------------------------'
print ' (0) - Import VTK OFF'
print ' (1) - Import VTK ON'
import_option = int(raw_input("\n enter option above: "))
if import_option == 1:
 folderName = raw_input("\n enter simulation folder name VTK import: ")
 numberStep = int(raw_input("\n enter number step VTK import: "))
print' ----------------------------------------------------------------------------\n'



print ' ----------------------------------------------------------------------------'
print ' (1) - ALE'
print ' (2) - Eulerian'
description_option = int(raw_input("\n enter motion description option above: "))
if description_option == 1:
 kLagrangian = float(raw_input("\n enter the control parameter of the Lagrangian: "))
 kLaplacian = float(raw_input("\n enter the control parameter of the Laplacian Smoothing: "))
 kVelocity = float(raw_input("\n enter the control parameter of the Velocity Smoothing: "))
 description_name = 'ALE'
else:
 description_name = 'Eulerian'
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' (0) - Boundary Moving OFF'
print ' (1) - Boundary Moving ON'
boundary_moving_option = int(raw_input("\n enter option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' (1) - Taylor Galerkin Scheme'
print ' (2) - Semi Lagrangian Scheme'
scheme_option = int(raw_input("\n Enter simulation scheme option above: "))
if scheme_option == 1:
 scheme_name = 'Taylor Galerkin'
elif scheme_option == 2:
 scheme_name = 'Semi Lagrangian'
print' ----------------------------------------------------------------------------\n'



print ' ----------------------------------------------------------------------------'
print ' (0) - Analytic Linear Element'
print ' (1) - Linear Element'
print ' (2) - Mini Element'
print ' (3) - Quadratic Element'
print ' (4) - Cubic Element'
polynomial_option = int(raw_input("\n Enter polynomial degree option above: "))
print' ----------------------------------------------------------------------------\n'


if simulation_option == 1:
 if polynomial_option == 0:
  gausspoints = 3

 else:
  print ' ----------------------------------------------------------------------------'
  print ' 3 Gauss Points'
  print ' 4 Gauss Points'
  print ' 6 Gauss Points'
  print ' 12 Gauss Points'
  gausspoints = int(raw_input("\n Enter Gauss Points Number option above: "))
  print' ----------------------------------------------------------------------------\n'


 
 print ' ----------------------------------------------------------------------------'
 nt = int(raw_input(" Enter number of time interations (nt): "))
 print' ----------------------------------------------------------------------------\n'
 
 
 print ' ----------------------------------------------------------------------------'
 folderResults = raw_input(" Enter folder name to save simulations: ")
 print' ----------------------------------------------------------------------------\n'

 print ' ----------------------------------------------------------------------------'
 observation = raw_input(" Digit observation: ")
 print' ----------------------------------------------------------------------------\n'


elif simulation_option == 0:
 gausspoints = 3
 nt = 20
 folderResults  = 'deletar'
 observation = 'debug'



print '\n ------------'
print ' IMPORT MESH:'
print ' ------------'


solution_start_time = time()

import_mesh_start_time = time()

# Linear and Mini Elements
if polynomial_option == 0 or polynomial_option == 1 or polynomial_option == 2:
 #mshFileName = 'linearHalfPoiseuille.msh'
 #mshFileName = 'linearStraightGeo.msh'
 #mshFileName = 'poiseuille.msh'
 #mshFileName = 'poiseuilleV2.msh'
 #mshFileName = 'mesh5.msh'
 mshFileName = 'CurvedGeoStrut.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 #Linear Element
 if polynomial_option == 0 or polynomial_option == 1:
  mesh = importMSH.Linear2D(pathMSHFile, mshFileName)

  numNodes               = mesh.numNodes
  numVerts               = mesh.numNodes
  numElements            = mesh.numElements
  x                      = mesh.x
  y                      = mesh.y
  IEN                    = mesh.IEN
  boundaryEdges          = mesh.boundaryEdges
  boundaryNodes          = mesh.boundaryNodes
  neighborsNodes         = mesh.neighborsNodes
  neighborsNodesALE      = mesh.neighborsNodesALE
  neighborsElements      = mesh.neighborsElements
  minLengthMesh          = mesh.minLengthMesh
  FreedomDegree          = mesh.FreedomDegree
  numPhysical            = mesh.numPhysical 

  Re = 1.0
  Sc = 1.0
  CFL = 0.5
  dt = float(CFL*minLengthMesh)
  #dt = 0.1   #SL 

 #Mini Element
 elif polynomial_option == 2:
  mesh = importMSH.Mini2D(pathMSHFile, mshFileName)

  numNodes                    = mesh.numNodes
  numVerts                    = mesh.numVerts
  numElements                 = mesh.numElements
  x                           = mesh.x
  y                           = mesh.y
  IEN                         = mesh.IEN
  boundaryEdges               = mesh.boundaryEdges
  boundaryNodes               = mesh.boundaryNodes
  neighborsNodes              = mesh.neighborsNodes
  neighborsNodesALE           = mesh.neighborsNodesALE
  neighborsNodesPressure      = mesh.neighborsNodesPressure
  neighborsElements           = mesh.neighborsElements
  minLengthMesh               = mesh.minLengthMesh
  velocityFreedomDegree       = mesh.velocityFreedomDegree
  pressureFreedomDegree       = mesh.pressureFreedomDegree
  numPhysical                 = mesh.numPhysical 
  Re = 54.5
  Sc = 10.0
  CFL = 0.5
  dt = float(CFL*minLengthMesh)
  #dt = 0.01   #linear result ok 


# Quad Element
elif polynomial_option == 3:
 #mshFileName = 'quadHalfPoiseuille.msh'
 #mshFileName = 'quadStraightGeo.msh'
 #mshFileName = 'quadCurvedGeoStrut.msh'
 mshFileName = 'quadRealGeoStrut.msh'
 #mshFileName = 'quadCurvedGeo.msh'
 #mshFileName = 'quad1.msh'
 #mshFileName = 'quad2.msh'
 #mshFileName = 'quad3.msh'
 #mshFileName = 'quad4.msh'
 #mshFileName = 'quad5.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Quad2D(pathMSHFile, mshFileName)

 numNodes               = mesh.numNodes
 numVerts               = mesh.numVerts
 numElements            = mesh.numElements
 x                      = mesh.x
 y                      = mesh.y
 IEN                    = mesh.IEN
 boundaryEdges          = mesh.boundaryEdges
 boundaryNodes          = mesh.boundaryNodes
 neighborsNodes         = mesh.neighborsNodes
 neighborsNodesALE      = mesh.neighborsNodesALE
 neighborsElements      = mesh.neighborsElements
 minLengthMesh          = mesh.minLengthMesh
 FreedomDegree          = mesh.FreedomDegree
 numPhysical            = mesh.numPhysical 

 Re = 1.0
 Sc = 1.0
 CFL = 1.0
 vMax = 4.0
 dt = float(CFL*minLengthMesh/vMax)
 #dt = 0.0005  




# Cubic Element
elif polynomial_option == 4:
 mshFileName = 'cubicStent_cubic.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Cubic2D(pathMSHFile, mshFileName, equation_number)
 mesh.coord()
 mesh.ien()


import_mesh_end_time = time()
import_mesh_time = import_mesh_end_time - import_mesh_start_time
print ' time duration: %.1f seconds \n' %import_mesh_time



print ' ---------'
print ' ASSEMBLY:'
print ' ---------'


assembly_start_time = time()

Kxxr, Kxyr, Kyxr, Kyyr, Kr, M2r, Mr, M, MrLump, Gx, Gy, Gxr, Gyr, M1, polynomial_order = assembly.AxiNS2D(simulation_option, polynomial_option, velocityFreedomDegree, pressureFreedomDegree, numNodes, numVerts, numElements, IEN, x, y, gausspoints)


#prof. gustavo reference
#A11 = (Mr/dt) + (1./Re)*(Kxxr + Kyyr + 2.0*M)
#A13 = -Gxr
#A22 = (Mr/dt) + (1./Re)*(Kxxr + Kyyr)
#A23 = -Gyr
#A31 = Gxr.transpose() + M1
#A32 = Gyr.transpose()
#
#A11 = sps.csr_matrix.tolil(A11)         
#A13 = sps.csr_matrix.tolil(A13)         
#A22 = sps.csr_matrix.tolil(A22)         
#A23 = sps.csr_matrix.tolil(A23)         
#A31 = sps.csr_matrix.tolil(A31)         
#A32 = sps.lil_matrix.tolil(A32)         
#A = sps.bmat([[A11  , None , A13],              
#              [None , A22  , A23],
#              [A31  , A32  , None]], format='lil')             
#

# -----

#calculated using r_avg in MEF integral - result ok (COBEM 2021)
A11 = (Mr/dt) + (1./Re)*(Kyyr + Kxxr)
A13 = -Gxr
A22 = (Mr/dt) + (1./Re)*(Kyyr + Kxxr)
A23 = -Gyr
A31 = Gxr.transpose()
A32 = Gyr.transpose()

A11 = sps.csr_matrix.tolil(A11)         
A13 = sps.csr_matrix.tolil(A13)         
A22 = sps.csr_matrix.tolil(A22)         
A23 = sps.csr_matrix.tolil(A23)         
A31 = sps.lil_matrix.tolil(A31)         
A32 = sps.lil_matrix.tolil(A32)         
A = sps.bmat([[A11  , None , A13],              
              [None , A22  , A23],
              [A31  , A32  , None]], format='lil')             


concentrationLHS = (np.copy(Mr)/dt) + (1.0/(Re*Sc))*np.copy(Kxxr) + (1.0/(Re*Sc))*np.copy(Kyyr)

#concentrationLHS = concentrationLHS.todense()
#concentrationLHS = np.array(concentrationLHS, dtype = float)


assembly_end_time = time()
assembly_time = assembly_end_time - assembly_start_time
print ' time duration: %.1f seconds \n' %assembly_time





print ' --------------------------------'
print ' INITIAL AND BOUNDARY CONDITIONS:'
print ' --------------------------------'

bc_apply_start_time = time()


# ------------------------ Boundaries Conditions ----------------------------------

# Linear and Mini Elements
if polynomial_option == 0 or polynomial_option == 1 or polynomial_option == 2:

 # Applying vx condition
 xVelocityBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
 xVelocityBC.xVelocityCondition(boundaryEdges,neighborsNodes)
 benchmark_problem = xVelocityBC.benchmark_problem

 # Applying vy condition
 yVelocityBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
 yVelocityBC.yVelocityCondition(boundaryEdges,neighborsNodes)
 
 # Applying pressure condition
 pressureBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
 pressureBC.pressureCondition(boundaryEdges,neighborsNodesPressure)

 # Applying concentration condition
 concentrationBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes, numVerts,x,y)
 concentrationBC.concentrationCondition(boundaryEdges,neighborsNodes)

 for mm in xVelocityBC.dirichletNodes:
   A[mm,:] = 0.0 
   A[mm,mm] = 1.0 

 for mm in yVelocityBC.dirichletNodes:
   A[mm + numNodes,:] = 0.0 
   A[mm + numNodes,mm + numNodes] = 1.0 

 for mm in pressureBC.dirichletNodes:
   A[mm + 2*numNodes,:] = 0.0 
   A[mm + 2*numNodes,mm + 2*numNodes] = 1.0 

 for mm in concentrationBC.dirichletNodes:
   concentrationLHS[mm,:] = 0.0 
   concentrationLHS[mm,mm] = 1.0 




# Quad Element
elif polynomial_option == 3:

 # Applying vx condition
 xVelocityLHS0 = sps.lil_matrix.copy(M)
 xVelocityBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
 xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
 benchmark_problem = xVelocityBC.benchmark_problem

 # Applying vr condition
 yVelocityLHS0 = sps.lil_matrix.copy(M)
 yVelocityBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
 yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
 streamFunctionBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
 streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 vorticityDirichletNodes = boundaryNodes

 # Applying concentration condition
 concentrationLHS0 = (sps.lil_matrix.copy(M)/dt) + (1.0/(Re*Sc))*sps.lil_matrix.copy(Kxx) + (1.0/(Re*Sc))*sps.lil_matrix.copy(Kyy)
 concentrationBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
 concentrationBC.concentrationCondition(boundaryEdges,concentrationLHS0,neighborsNodes)
# ---------------------------------------------------------------------------------




# -------------------------- Import VTK File - Initial Condition ------------------------------------
if import_option == 0:
 import_option = 'OFF'
 
 vx = np.copy(xVelocityBC.aux1BC)
 vy = np.copy(yVelocityBC.aux1BC)
 p = np.copy(pressureBC.aux1BC)
 c = np.copy(concentrationBC.aux1BC)
 sol = np.concatenate((vx, vy, p), axis=0)

elif import_option == 1:
 import_option = 'ON'
 
 vx, vy, c, p, p, benchmark_problem = importVTK.vtkFile("/home/marquesleandro/NSAxi/results/" + folderName + "/" + folderName + str(numberStep) + ".vtk", numNodes, numVerts, numElements, IEN, polynomial_option)
 sol = np.concatenate((vx, vy, p), axis=0)
#----------------------------------------------------------------------------------


bc_apply_end_time = time()
bc_apply_time = bc_apply_end_time - bc_apply_start_time 
print ' time duration: %.1f seconds \n' %bc_apply_time







print ' -----------------------------'
print ' PARAMETERS OF THE SIMULATION:'
print ' -----------------------------'

print ' Benchmark Problem: %s' %benchmark_problem
print ' Motion Description: %s' %str(description_name)
print ' Scheme: %s' %str(scheme_name)
print ' Element Type: %s' %str(polynomial_order)
print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
print ' Mesh: %s' %mshFileName
print ' Number of nodes: %s' %numNodes
print ' Number of elements: %s' %numElements
print ' Smallest edge length: %f' %minLengthMesh
print ' Time step: %s' %dt
print ' Import VTK: %s' %import_option
print ' Number of time iteration: %s' %nt
print ' Reynolds number: %s' %Re
print ' Schmidt number: %s' %Sc
print ""


print ' ----------------------------'
print ' SOLVE THE LINEARS EQUATIONS:'
print ' ---------------------------- \n'

print ' Saving simulation in %s \n' %folderResults



os.chdir(initial_path)



# ------------------------ Export VTK File ---------------------------------------
# Linear and Mini Elements
if polynomial_option == 0 or polynomial_option == 1 or polynomial_option == 2:   
 save = exportVTK.Linear2D(x,y,IEN,numVerts,numElements,c,p,p,vx,vy)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))

# Quad Element
elif polynomial_option == 3:   
 save = exportVTK.Quad2D(x,y,IEN,numNodes,numElements,w,psi,c,vx,vy)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))
# ---------------------------------------------------------------------------------



x_old = np.zeros([numNodes,1], dtype = float)
y_old = np.zeros([numNodes,1], dtype = float)
vx_old = np.zeros([numNodes,1], dtype = float)
vy_old = np.zeros([numNodes,1], dtype = float)
end_type = 0
for t in tqdm(range(1, nt)):
 numIteration = t

 try:
  print ""
  print '''
                 COPYRIGHT                    
   ===========================================
   Simulator: %s
   created by Leandro Marques at 02/2021
   e-mail: marquesleandro67@gmail.com
   COPPE/Departament of Mechanical Engineering
   Federal University of the Rio de Janeiro
   ===========================================
  ''' %sys.argv[0]
 
 
 
  print ' -----------------------------'
  print ' PARAMETERS OF THE SIMULATION:'
  print ' -----------------------------'
 
  print ' Benchmark Problem: %s' %benchmark_problem
  print ' Motion Description: %s' %str(description_name)
  print ' Scheme: %s' %str(scheme_name)
  print ' Element Type: %s' %str(polynomial_order)
  print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
  print ' Mesh: %s' %mshFileName
  print ' Number of nodes: %s' %numNodes
  print ' Number of elements: %s' %numElements
  print ' Smallest edge length: %f' %minLengthMesh
  print ' Time step: %s' %dt
  print ' Import VTK: %s' %import_option
  print ' Number of time iteration: %s' %numIteration
  print ' Reynolds number: %s' %Re
  print ' Schmidt number: %s' %Sc


 
  # ------------------------- Boundary Moving ------------------------------------

  if boundary_moving_option == 1:
   y_boundary = np.zeros([numNodes,1], dtype = float)
 
   for i in range(0,len(boundaryEdges)):
    line = boundaryEdges[i][0]
    v1 = boundaryEdges[i][1] - 1
    v2 = boundaryEdges[i][2] - 1
  
    # oscillatory parameters ok
    # center point x=5.3
    if line == 1 or line == 5:
     y_boundary[v1] = 2.0*0.01*np.sin((2.0*np.pi/7.0)*x[v1])*np.cos((2.0*np.pi/32.0)*t)
     y_boundary[v2] = 2.0*0.01*np.sin((2.0*np.pi/7.0)*x[v2])*np.cos((2.0*np.pi/32.0)*t)
  
   y = y + y_boundary
  # --------------------------------------------------------------------------------
 

 
  # ------------------------- ALE Scheme --------------------------------------------

  ALE_start_time = time()

  if description_option == 1:
   xmeshALE_dif = np.linalg.norm(x-x_old)
   ymeshALE_dif = np.linalg.norm(y-y_old)
   if not xmeshALE_dif < 5e-3 and not ymeshALE_dif < 5e-3:
    x_old = np.copy(x)
    y_old = np.copy(y)
   
    print ""
    print ' ------------'
    print ' MESH UPDATE:'
    print ' ------------'
   
   
    vxLaplacianSmooth, vyLaplacianSmooth = ALE.Laplacian_smoothing(neighborsNodesALE, numNodes, numVerts, x, y, dt)
    #vxLaplacianSmooth, vyLaplacianSmooth = ALE.MINILaplacian_smoothing(neighborsNodesALE, numNodes, numVerts, numElements, IEN, x, y, dt)
    #vxLaplacianSmooth, vyLaplacianSmooth = ALE.Laplacian_smoothing_avg(neighborsNodesALE, numNodes, numVerts, x, y, dt)
    vxVelocitySmooth,  vyVelocitySmooth  = ALE.Velocity_smoothing(neighborsNodes, numNodes, vx, vy)
  
    vxMesh = kLagrangian*vx + kLaplacian*vxLaplacianSmooth + kVelocity*vxVelocitySmooth
    vyMesh = kLagrangian*vy + kLaplacian*vyLaplacianSmooth + kVelocity*vyVelocitySmooth
   
   
    for i in boundaryNodes:
     node = i-1 
     vxMesh[node] = 0.0
     vyMesh[node] = 0.0
   
    x = x + vxMesh*dt
    y = y + vyMesh*dt

    x = np.asarray(x) 
    y = np.asarray(y)

    for e in range(0,numElements):
     v1 = IEN[e][0]
     v2 = IEN[e][1]
     v3 = IEN[e][2]
     v4 = IEN[e][3]
    
     vxMesh[v4] = (vxMesh[v1] + vxMesh[v2] + vxMesh[v3])/3.0
     vyMesh[v4] = (vyMesh[v1] + vyMesh[v2] + vyMesh[v3])/3.0
     x[v4] = (x[v1] + x[v2] + x[v3])/3.0
     y[v4] = (y[v1] + y[v2] + y[v3])/3.0
    

     
    vxALE = vx - vxMesh
    vyALE = vy - vyMesh
   
    ALE_end_time = time()
    ALE_time_solver = ALE_end_time - ALE_start_time
    print ' time duration: %.1f seconds' %ALE_time_solver
    # ---------------------------------------------------------------------------------
   
  
  
  
    # ------------------------- Assembly --------------------------------------------
    print ""
    print ' ---------'
    print ' ASSEMBLY:'
    print ' ---------'
  
    Kxxr, Kxyr, Kyxr, Kyyr, Kr, M2r, Mr, M, MrLump, Gx, Gy, Gxr, Gyr, M1, polynomial_order = assembly.AxiNS2D(simulation_option, polynomial_option, velocityFreedomDegree, pressureFreedomDegree, numNodes, numVerts, numElements, IEN, x, y, gausspoints)
    
    
    #prof. gustavo reference
    #A11 = (Mr/dt) + (1./Re)*(Kxxr + Kyyr + 2.0*M)
    #A13 = -Gxr
    #A22 = (Mr/dt) + (1./Re)*(Kxxr + Kyyr)
    #A23 = -Gyr
    #A31 = Gxr.transpose() + M1
    #A32 = Gyr.transpose()
    #
    #A11 = sps.csr_matrix.tolil(A11)         
    #A13 = sps.csr_matrix.tolil(A13)         
    #A22 = sps.csr_matrix.tolil(A22)         
    #A23 = sps.csr_matrix.tolil(A23)         
    #A31 = sps.csr_matrix.tolil(A31)         
    #A32 = sps.lil_matrix.tolil(A32)         
    #A = sps.bmat([[A11  , None , A13],              
    #              [None , A22  , A23],
    #              [A31  , A32  , None]], format='lil')             
    #
    
    # -----
    
    #calculated using r_avg in MEF integral - result ok (COBEM 2021)
    A11 = (Mr/dt) + (1./Re)*(Kyyr + Kxxr)
    A13 = -Gxr
    A22 = (Mr/dt) + (1./Re)*(Kyyr + Kxxr)
    A23 = -Gyr
    A31 = Gxr.transpose()
    A32 = Gyr.transpose()
    
    A11 = sps.csr_matrix.tolil(A11)         
    A13 = sps.csr_matrix.tolil(A13)         
    A22 = sps.csr_matrix.tolil(A22)         
    A23 = sps.csr_matrix.tolil(A23)         
    A31 = sps.lil_matrix.tolil(A31)         
    A32 = sps.lil_matrix.tolil(A32)         
    A = sps.bmat([[A11  , None , A13],              
                  [None , A22  , A23],
                  [A31  , A32  , None]], format='lil')             
    
    
    concentrationLHS = (np.copy(Mr)/dt) + (1.0/(Re*Sc))*np.copy(Kxxr) + (1.0/(Re*Sc))*np.copy(Kyyr)
    
    #concentrationLHS = concentrationLHS.todense()
    #concentrationLHS = np.array(concentrationLHS, dtype = float)

    # --------------------------------------------------------------------------------
  
  
  
  
    # ------------------------ Boundaries Conditions ----------------------------------
    print ""
    print ' --------------------------------'
    print ' INITIAL AND BOUNDARY CONDITIONS:'
    print ' --------------------------------'
 
    start_BC_time = time()

    # Linear and Mini Elements
    if polynomial_option == 0 or polynomial_option == 1 or polynomial_option == 2:
    
     # Applying vx condition
     xVelocityBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
     xVelocityBC.xVelocityCondition(boundaryEdges,neighborsNodes)
     benchmark_problem = xVelocityBC.benchmark_problem
    
     # Applying vy condition
     yVelocityBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
     yVelocityBC.yVelocityCondition(boundaryEdges,neighborsNodes)
     
     # Applying pressure condition
     pressureBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes,numVerts,x,y)
     pressureBC.pressureCondition(boundaryEdges,neighborsNodesPressure)
    
     # Applying concentration condition
     concentrationBC = benchmarkProblems.AxiNS2DStent(numPhysical,numNodes, numVerts,x,y)
     concentrationBC.concentrationCondition(boundaryEdges,neighborsNodes)
    
     for mm in xVelocityBC.dirichletNodes:
       A[mm,:] = 0.0 
       A[mm,mm] = 1.0 
    
     for mm in yVelocityBC.dirichletNodes:
       A[mm + numNodes,:] = 0.0 
       A[mm + numNodes,mm + numNodes] = 1.0 
    
     for mm in pressureBC.dirichletNodes:
       A[mm + 2*numNodes,:] = 0.0 
       A[mm + 2*numNodes,mm + 2*numNodes] = 1.0 
    
     for mm in concentrationBC.dirichletNodes:
       concentrationLHS[mm,:] = 0.0 
       concentrationLHS[mm,mm] = 1.0 
    
    
    
    
    # Quad Element
    elif polynomial_option == 3:
    
     # Applying vx condition
     xVelocityLHS0 = sps.lil_matrix.copy(M)
     xVelocityBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
     xVelocityBC.xVelocityCondition(boundaryEdges,xVelocityLHS0,neighborsNodes)
     benchmark_problem = xVelocityBC.benchmark_problem
    
     # Applying vr condition
     yVelocityLHS0 = sps.lil_matrix.copy(M)
     yVelocityBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
     yVelocityBC.yVelocityCondition(boundaryEdges,yVelocityLHS0,neighborsNodes)
     
     # Applying psi condition
     streamFunctionLHS0 = sps.lil_matrix.copy(Kxx) + sps.lil_matrix.copy(Kyy)
     streamFunctionBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
     streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)
    
     # Applying vorticity condition
     vorticityDirichletNodes = boundaryNodes
    
     # Applying concentration condition
     concentrationLHS0 = (sps.lil_matrix.copy(M)/dt) + (1.0/(Re*Sc))*sps.lil_matrix.copy(Kxx) + (1.0/(Re*Sc))*sps.lil_matrix.copy(Kyy)
     concentrationBC = benchmarkProblems.quadStent(numPhysical,numNodes,x,y)
     concentrationBC.concentrationCondition(boundaryEdges,concentrationLHS0,neighborsNodes)


    end_BC_time = time()
    bc_apply_time = end_BC_time - start_BC_time
    print ' time duration: %.1f seconds \n' %bc_apply_time

  elif description_option == 2:
   start_time = time()
   vxMesh = 0.0*vx
   vyMesh = 0.0*vy

   vxALE = vx - vxMesh
   vyALE = vy - vyMesh

   ALE_end_time = time()
   ALE_time_solver = ALE_end_time - ALE_start_time
   print ' time duration: %.1f seconds' %ALE_time_solver
 
  # ---------------------------------------------------------------------------------
   
   
 



  # ------------------------ semi-Lagrangian Method --------------------------------
  if scheme_option == 2:
   print ""
   print ' -----------------------'
   print ' SEMI-LAGRANGIAN METHOD:'
   print ' -----------------------'
   start_SL_time = time()

   # Linear Element   
   if polynomial_option == 0 or polynomial_option == 1:
    vx_d, vy_d = semiLagrangian.Linear2D(numNodes, neighborsElements, IEN, x, y, vxALE, vyALE, dt, vx, vy)

   # Mini Element   
   elif polynomial_option == 2:
    vx_d, vy_d, c_d = semiLagrangian.Mini2D(numNodes, neighborsElements, IEN, x, y, vxALE, vyALE, dt, vx, vy, c)
 
   # Quad Element   
   elif polynomial_option == 3:
    w_d, c_d = semiLagrangian.Quad2D(numNodes, neighborsElements, IEN, x, y, vxALE, vyALE, dt, w, c)
 
   end_SL_time = time()
   SL_time = end_SL_time - start_SL_time
   print ' time duration: %.1f seconds' %SL_time
  #----------------------------------------------------------------------------------






  # ------------------------ SOLVE LINEAR EQUATIONS ----------------------------------
  print ""
  print ' ----------------------------'
  print ' SOLVE THE LINEARS EQUATIONS:'
  print ' ----------------------------'

  #---------- Step 1 - Solve the continuity and momentum equation ----------------------
  start_solver_time = time()
  start_ns_time = time()

#  AA = np.copy(M2r.todense()/dt)
#  b = np.dot(AA,np.concatenate((vx_d,vy_d),axis=0))
#  #b = np.dot(M/dt,np.concatenate((vx,vy),axis=0))  #stokes
#  bp = np.zeros([numVerts,1], dtype = float)
#  b = np.concatenate((b,bp),axis=0)
#
#  for i in xVelocityBC.dirichletNodes:
#   b[i] = xVelocityBC.aux1BC[i]
#
#  for i in yVelocityBC.dirichletNodes:
#   b[i + numNodes] = yVelocityBC.aux1BC[i]
#
#  for i in pressureBC.dirichletNodes:
#   b[i + 2*numNodes] = pressureBC.aux1BC[i]
#
#  sol = scipy.sparse.linalg.spsolve(A.tocsr(),b)
#  #sol = scipy.sparse.linalg.cg(A,b, sol, maxiter=1.0e+03, tol=1.0e-04)
#  #sol = sol[0].reshape((len(sol[0]),1))
#
#  vx = np.reshape(np.array(sol[0:numNodes]),(len(vx),1))
#  vy = np.reshape(np.array(sol[numNodes:2*numNodes]),(len(vy),1))
#  p  = np.reshape(np.array(sol[2*numNodes:]),(len(p),1))
#
#  end_ns_time = time()
#  ns_time = end_ns_time - start_ns_time
#  print ' Solver Continuity and Momentum Equation: %.1f seconds' %ns_time
#  #----------------------------------------------------------------------------------
# 
# 
#
#  #---------- Step 2 - Solve the specie transport equation ----------------------
#  start_concentrationsolver_time = time()
#
#  c_old = np.copy(c)
#  # Taylor Galerkin Scheme
#  if scheme_option == 1:
#   A = np.copy(M)/dt 
#   concentrationRHS = sps.lil_matrix.dot(A,c) - np.multiply(vx,sps.lil_matrix.dot(Gx,c))\
#         - np.multiply(vy,sps.lil_matrix.dot(Gy,c))\
#         - (dt/2.0)*np.multiply(vx,(np.multiply(vx,sps.lil_matrix.dot(Kxx,c)) + np.multiply(vy,sps.lil_matrix.dot(Kyx,c))))\
#         - (dt/2.0)*np.multiply(vy,(np.multiply(vx,sps.lil_matrix.dot(Kxy,c)) + np.multiply(vy,sps.lil_matrix.dot(Kyy,c))))
#   concentrationRHS = np.multiply(concentrationRHS,concentrationBC.aux2BC)
#   concentrationRHS = concentrationRHS + concentrationBC.dirichletVector
#   c = scipy.sparse.linalg.cg(concentrationBC.LHS,concentrationRHS,c, maxiter=1.0e+05, tol=1.0e-05)
#   c = c[0].reshape((len(c[0]),1))
# 
# 
# 
#  # Semi-Lagrangian Scheme
#  elif scheme_option == 2:
#   concentrationRHS = np.dot(Mr.todense()/dt,c_d)
# 
#   for i in concentrationBC.dirichletNodes:
#    concentrationRHS[i] = concentrationBC.aux1BC[i]
#
#   c = scipy.sparse.linalg.spsolve(concentrationLHS.tocsr(),concentrationRHS)
#   #c = scipy.sparse.linalg.cg(concentrationLHS,concentrationRHS, c, maxiter=1.0e+05, tol=1.0e-05)
#   #c = c[0].reshape((len(c[0]),1))
#
#
#
#
#  end_concentrationsolver_time = time()
#  concentrationsolver_time = end_concentrationsolver_time - start_concentrationsolver_time
#  print ' Concentration Solver: %.1f seconds' %concentrationsolver_time
#
#
#
  end_solver_time = time()
  solver_time = end_solver_time - start_solver_time
  #----------------------------------------------------------------------------------
 

 
 
 
 

  # ------------------------ Export VTK File ---------------------------------------
  print ""
  print ' ----------------'
  print ' EXPORT VTK FILE:'
  print ' ----------------'
  print ' Saving simulation in %s' %folderResults
  start_exportVTK_time = time()

  # Linear and Mini Elements
  if polynomial_option == 0 or polynomial_option == 1 or polynomial_option == 2:   
   save = exportVTK.Linear2D(x,y,IEN,numVerts,numElements,c,p,p,vx,vy)
   save.create_dir(folderResults)
   save.saveVTK(folderResults + str(t))
 
  # Quad Element
  elif polynomial_option == 3:   
   save = exportVTK.Quad2D(x,y,IEN,numNodes,numElements,w,psi,c,vx,vy)
   save.create_dir(folderResults)
   save.saveVTK(folderResults + str(t))

  end_exportVTK_time = time()
  exportVTK_time = end_exportVTK_time - start_exportVTK_time
  print ' time duration: %.1f seconds' %exportVTK_time
  #----------------------------------------------------------------------------------
 




 
 
  # ---------------------------------------------------------------------------------
  print ""
  print ' -------'
  print ' CHECKS:'
  print ' -------'
  start_checks_time = time()
 
  # CHECK STEADY STATE
  #if np.all(vx == vx_old) and np.all(vy == vy_old):
  # end_type = 1
  # break
 
  # CHECK CONVERGENCE OF THE SOLUTION
  if np.linalg.norm(vx) > 10e2 or np.linalg.norm(vy) > 10e2:
   end_type = 2
   break

  end_checks_time = time()
  checks_time = end_checks_time - start_checks_time
  print ' time duration: %.1f seconds' %checks_time
  # ---------------------------------------------------------------------------------
 
  print "" 
  print "" 
  print " ---------------------------------------------------------------------------------"
  


 except KeyboardInterrupt:
  end_type = 3
  break 

 # used to calculate code performance
 if t == 1: 
  end_first_loop_time = time()
  first_loop_time = end_first_loop_time - solution_start_time


end_time = time()
solution_time = end_time - solution_start_time


print ""
print ' ----------------'
print ' SAVING RELATORY:'
print ' ----------------'

if end_type == 0:
 print ' END SIMULATION. NOT STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 1:
 print ' END SIMULATION. STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 2:
 print ' END SIMULATION. ERROR CONVERGENCE RESULT'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 3:
 print ' END SIMULATION. FORCED INTERRUPTION'
 print ' Relatory saved in %s' %folderResults
 print ""




# -------------------------------- Export Relatory ---------------------------------------
relatory.export(save.path, folderResults, sys.argv[0], benchmark_problem, description_name, scheme_name, mshFileName, numNodes, numElements, minLengthMesh, dt, numIteration, Re, Sc, import_mesh_time, assembly_time, bc_apply_time, ALE_time_solver, SL_time, solver_time, first_loop_time, solution_time, polynomial_order, gausspoints, observation)
# ----------------------------------------------------------------------------------------



