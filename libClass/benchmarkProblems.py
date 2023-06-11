# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to compute boundary condition


import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

# OBS: para vetor devemos unir eixo x e y no mesmo vetor, logo usar np.row_stack([dirichlet_pts[1],dirichlet_pts[2]])


class linearStent:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'linear Stent'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v1]**2)
    _self.aux1BC[v2] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v2]**2)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Symmetric axis
   elif line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Symmetric axis (Bottom Line)
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Noslip (Top Line)
   # Ref: Batchelor 1967 pag. 78 eq. 2.2.12
   # As psi_bottom is zero, so psi_top is:
   elif line == 1 or line == 5:
    _self.aux1BC[v1] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v2] = (_self.maxVx*(2.0/3.0))*(_self.L)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 

 def concentrationCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Stent
   if line == 5:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 





class quadStent:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'quad Stent'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Noslip 
   if line == 1 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v1]**2)
    _self.aux1BC[v2] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v2]**2)
    _self.aux1BC[v3] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v3]**2)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Noslip 
   if line == 1 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Symmetric axis
   elif line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Symmetric axis (Bottom Line)
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Noslip (Top Line)
   # Ref: Batchelor 1967 pag. 78 eq. 2.2.12
   # As psi_bottom is zero, so psi_top is:
   elif line == 1 or line == 5:
    _self.aux1BC[v1] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v2] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v3] = (_self.maxVx*(2.0/3.0))*(_self.L)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 

 def concentrationCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Stent
   if line == 5:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0
    _self.aux1BC[v3] = 1.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



class Convection2D:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying c condition
 # condition_concentration = bc_apply.Convection(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_concentration.neumann_condition(mesh.neumann_edges[1])
 # condition_concentration.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_concentration.gaussian_elimination(LHS_c0,mesh.neighbors_nodes)
 # condition_concentration.initial_condition()
 # c = np.copy(condition_concentration.c)
 # vx = np.copy(condition_concentration.vx)
 # vy = np.copy(condition_concentration.vy)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _nphysical, _npoints, _x, _y):
  _self.nphysical = _nphysical
  _self.npoints = _npoints
  _self.x = _x
  _self.y = _y
  _self.bc = np.zeros([_self.nphysical,1], dtype = float) 
  _self.benchmark_problem = 'Convection 2D'

  # Velocity c condition
  _self.bc[0][0] = 0.0
  _self.bc[1][0] = 0.0
  _self.bc[2][0] = 0.0
  _self.bc[3][0] = 0.0


 def neumann_condition(_self, _neumann_edges):
  _self.bc_neumann = np.zeros([_self.npoints,1], dtype = float) 
  _self.neumann_edges = _neumann_edges 
 
  for i in range(0, len(_self.neumann_edges)):
   line = _self.neumann_edges[i][0] - 1
   v1 = _self.neumann_edges[i][1] - 1
   v2 = _self.neumann_edges[i][2] - 1

   x = _self.x[v1] - _self.x[v2]
   y = _self.y[v1] - _self.y[v2]
   length = np.sqrt(x**2 + y**2)
  
   _self.bc_neumann[v1] += (_self.bc[line]*length) / 2. 
   _self.bc_neumann[v2] += (_self.bc[line]*length) / 2. 


 def dirichlet_condition(_self, _dirichlet_pts):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.dirichlet_pts = _dirichlet_pts
 

  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1

   _self.bc_1[v1] = _self.bc[line]
   _self.bc_1[v2] = _self.bc[line]

   _self.bc_neumann[v1] = 0.0 #Dirichlet condition is preferential
   _self.bc_neumann[v2] = 0.0 #Dirichlet condition is preferential

   _self.ibc.append(v1)
   _self.ibc.append(v2)
   
  _self.ibc = np.unique(_self.ibc)


 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.neighbors_nodes = _neighbors_nodes

  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 

 def initial_condition(_self):
  _self.c = np.copy(_self.bc_dirichlet)
  _self.vx = np.zeros([_self.npoints,1], dtype = float)
  _self.vy = np.zeros([_self.npoints,1], dtype = float)

  for i in range(0,_self.npoints):
   _self.vx[i] = _self.y[i]/3.0
   _self.vy[i] = -_self.x[i]/3.0

  a = 0.0
  b = -1.5
  r = 1.0

  for i in range(0, _self.npoints):
   x = _self.x[i] - a
   y = _self.y[i] - b
   lenght = np.sqrt(x**2 + y**2)

   if lenght < r:
    _self.c[i] = r**2 - (_self.x[i] - a)**2 - (_self.y[i] - b)**2


class Wave2D:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying c condition
 # condition_concentration = bc_apply.Convection(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_concentration.neumann_condition(mesh.neumann_edges[1])
 # condition_concentration.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_concentration.gaussian_elimination(LHS_c0,mesh.neighbors_nodes)
 # condition_concentration.initial_condition()
 # c = np.copy(condition_concentration.c)
 # vx = np.copy(condition_concentration.vx)
 # vy = np.copy(condition_concentration.vy)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _nphysical, _npoints, _x, _y):
  _self.nphysical = _nphysical
  _self.npoints = _npoints
  _self.x = _x
  _self.y = _y
  _self.bc = np.zeros([_self.nphysical,1], dtype = float) 
  _self.benchmark_problem = 'Wave 2D'

  # Velocity c condition
  _self.bc[0][0] = 0.0
  _self.bc[1][0] = 0.0
  _self.bc[2][0] = 0.0
  _self.bc[3][0] = 0.0


 def neumann_condition(_self, _neumann_edges):
  _self.bc_neumann = np.zeros([_self.npoints,1], dtype = float) 
  _self.neumann_edges = _neumann_edges 
 
  for i in range(0, len(_self.neumann_edges)):
   line = _self.neumann_edges[i][0] - 1
   v1 = _self.neumann_edges[i][1] - 1
   v2 = _self.neumann_edges[i][2] - 1

   x = _self.x[v1] - _self.x[v2]
   y = _self.y[v1] - _self.y[v2]
   length = np.sqrt(x**2 + y**2)
  
   _self.bc_neumann[v1] += (_self.bc[line]*length) / 2. 
   _self.bc_neumann[v2] += (_self.bc[line]*length) / 2. 


 def dirichlet_condition(_self, _dirichlet_pts):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.dirichlet_pts = _dirichlet_pts
 

  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1

   _self.bc_1[v1] = _self.bc[line]
   _self.bc_1[v2] = _self.bc[line]

   _self.bc_neumann[v1] = 0.0 #Dirichlet condition is preferential
   _self.bc_neumann[v2] = 0.0 #Dirichlet condition is preferential

   _self.ibc.append(v1)
   _self.ibc.append(v2)
   
  _self.ibc = np.unique(_self.ibc)


 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.neighbors_nodes = _neighbors_nodes

  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 

 def initial_condition(_self):
  _self.c = np.copy(_self.bc_dirichlet)
  _self.vx = np.zeros([_self.npoints,1], dtype = float)
  _self.vy = np.zeros([_self.npoints,1], dtype = float)

  for i in range(0,_self.npoints):
   _self.vx[i] = 0.0
   _self.vy[i] = 0.0

  a = 0.0
  b = 0.0
  r = 1.0

  for i in range(0, _self.npoints):
   x = _self.x[i] - a
   y = _self.y[i] - b
   lenght = np.sqrt(x**2 + y**2)

   if lenght < r:
    _self.c[i] = r**2 - (_self.x[i] - a)**2 - (_self.y[i] - b)**2







class linearHalfPoiseuille:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'linear Half Poiseuille'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v1]**2)
    _self.aux1BC[v2] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v2]**2)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Symmetric axis
   elif line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Symmetric axis (Bottom Line)
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Noslip (Top Line)
   # Ref: Batchelor 1967 pag. 78 eq. 2.2.12
   # As psi_bottom is zero, so psi_top is:
   elif line == 1:
    _self.aux1BC[v1] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v2] = (_self.maxVx*(2.0/3.0))*(_self.L)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 

 def concentrationCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Stent
   if line == 1:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 





class quadHalfPoiseuille:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'quad Half Poiseuille'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Noslip 
   if line == 1:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v1]**2)
    _self.aux1BC[v2] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v2]**2)
    _self.aux1BC[v3] = (_self.maxVx/(_self.L**2))*(_self.L**2 - _self.y[v3]**2)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Noslip 
   if line == 1:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Symmetric axis
   elif line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Symmetric axis (Bottom Line)
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
    _self.aux1BC[v3] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

   # Noslip (Top Line)
   # Ref: Batchelor 1967 pag. 78 eq. 2.2.12
   # As psi_bottom is zero, so psi_top is:
   elif line == 1:
    _self.aux1BC[v1] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v2] = (_self.maxVx*(2.0/3.0))*(_self.L)
    _self.aux1BC[v3] = (_self.maxVx*(2.0/3.0))*(_self.L)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 

 def concentrationCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1
   v3 = _self.boundaryEdges[i][3] - 1

   # Stent
   if line == 1:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0
    _self.aux1BC[v3] = 1.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)
    _self.dirichletNodes.append(v3)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



class NS2DPoiseuille:

 def __init__(_self, _numPhysical, _numNodes, _numVerts, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.numVerts = _numVerts
  _self.x = _x
  _self.y = _y
  _self.maxVx = 3.0/2.0
  _self.L = 1.0
  _self.benchmark_problem = 'NS2D Poiseuille'


 def xVelocityCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = (4.0*_self.maxVx*_self.y[v1])*(_self.L - _self.y[v1])
    _self.aux1BC[v2] = (4.0*_self.maxVx*_self.y[v2])*(_self.L - _self.y[v2])

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def yVelocityCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def pressureCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numVerts,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes


 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Outflow
   if line == 2:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def concentrationCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Stent
   if line == 4:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


class AxiNS2DStent:

 def __init__(_self, _numPhysical, _numNodes, _numVerts, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.numVerts = _numVerts
  _self.x = _x
  _self.y = _y
  _self.maxVx = 2.0
  _self.L = 1.0
  _self.benchmark_problem = 'NS2D Stent'


 def xVelocityCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = _self.maxVx*(1.0 - (_self.y[v1]/_self.L)**2)
    _self.aux1BC[v2] = _self.maxVx*(1.0 - (_self.y[v2]/_self.L)**2)

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def yVelocityCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 4 or line == 5:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Inflow
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def pressureCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numVerts,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes


 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Outflow
   if line == 2:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)



 def concentrationCondition(_self, _boundaryEdges, _neighborsNodes):
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Stent
   if line == 5:
    _self.aux1BC[v1] = 1.0
    _self.aux1BC[v2] = 1.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)






class QuadPoiseuille:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Poiseuille(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _nphysical, _npoints, _x, _y):
  _self.nphysical = _nphysical
  _self.npoints = _npoints
  _self.x = _x
  _self.y = _y
  _self.bc = np.zeros([_self.nphysical,1], dtype = float) 
  _self.benchmark_problem = 'Quad Poiseuille'

  # Velocity vx condition
  _self.bc[0][0] = 0.0
  _self.bc[1][0] = 0.0
  _self.bc[2][0] = 1.0
  _self.bc[3][0] = 0.0

  # Velocity vy condition
  _self.bc[4][0] = 0.0
  _self.bc[5][0] = 0.0
  _self.bc[6][0] = 0.0
  _self.bc[7][0] = 0.0


 def neumann_condition(_self, _neumann_edges):
  _self.bc_neumann = np.zeros([_self.npoints,1], dtype = float) 
  _self.neumann_edges = _neumann_edges 
 
  for i in range(0, len(_self.neumann_edges)):
   line = _self.neumann_edges[i][0] - 1
   v1 = _self.neumann_edges[i][1] - 1
   v2 = _self.neumann_edges[i][2] - 1
   v3 = _self.neumann_edges[i][3] - 1

   x1 = _self.x[v1] - _self.x[v3]
   y1 = _self.y[v1] - _self.y[v3]
   length1 = np.sqrt(x1**2 + y1**2)

   x2 = _self.x[v3] - _self.x[v2]
   y2 = _self.y[v3] - _self.y[v2]
   length2 = np.sqrt(x2**2 + y2**2)
 

   _self.bc_neumann[v1] += (_self.bc[line]*length1) / 2.0 
   _self.bc_neumann[v2] += (_self.bc[line]*length2) / 2.0 
   _self.bc_neumann[v3] += ((_self.bc[line]*length1) / 2.0) + ((_self.bc[line]*length2) / 2.0) 


 def dirichlet_condition(_self, _dirichlet_pts):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.dirichlet_pts = _dirichlet_pts
 

  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1
   v3 = _self.dirichlet_pts[i][3] - 1

   _self.bc_1[v1] = _self.bc[line]
   _self.bc_1[v2] = _self.bc[line]
   _self.bc_1[v3] = _self.bc[line]

   _self.bc_neumann[v1] = 0.0 #Dirichlet condition is preferential
   _self.bc_neumann[v2] = 0.0 #Dirichlet condition is preferential
   _self.bc_neumann[v3] = 0.0 #Dirichlet condition is preferential

   _self.ibc.append(v1)
   _self.ibc.append(v2)
   _self.ibc.append(v3)
   
  _self.ibc = np.unique(_self.ibc)


 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.neighbors_nodes = _neighbors_nodes

  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 


 def streamfunction_condition(_self, _dirichlet_pts, _LHS0, _neighbors_nodes):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.dirichlet_pts = _dirichlet_pts
  _self.neighbors_nodes = _neighbors_nodes

  # Dirichlet condition
  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1
   v2 = _self.dirichlet_pts[i][2] - 1
   v3 = _self.dirichlet_pts[i][3] - 1

   if line == 8:
    _self.bc_1[v1] = 0.0
    _self.bc_1[v2] = 0.0
    _self.bc_1[v3] = 0.0
 
    _self.ibc.append(v1)
    _self.ibc.append(v2)
    _self.ibc.append(v3)

   elif line == 11:
    _self.bc_1[v1] = 1.0
    _self.bc_1[v2] = 1.0
    _self.bc_1[v3] = 1.0

    _self.ibc.append(v1)
    _self.ibc.append(v2)
    _self.ibc.append(v3)

   elif line == 10:
    _self.bc_1[v1] = _self.y[v1]
    _self.bc_1[v2] = _self.y[v2]
    _self.bc_1[v3] = _self.y[v3]

    _self.ibc.append(v1)
    _self.ibc.append(v2)
    _self.ibc.append(v3)

  _self.ibc = np.unique(_self.ibc)


  # Gaussian elimination for psi
  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 





 
class Convection1D:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying c condition
 # condition_concentration = bc_apply.Convection(mesh.nphysical,mesh.npoints,mesh.x)
 # condition_concentration.neumann_condition(mesh.neumann_pts[1])
 # condition_concentration.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_concentration.gaussian_elimination(LHS_c0,mesh.neighbors_nodes)
 # condition_concentration.initial_condition()
 # c = np.copy(condition_concentration.c)
 # vx = np.copy(condition_concentration.vx)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _nphysical, _npoints, _x):
  _self.nphysical = _nphysical
  _self.npoints = _npoints
  _self.x = _x
  _self.bc = np.zeros([_self.nphysical,1], dtype = float) 
  _self.benchmark_problem = 'Convection 1D'

  # Velocity c condition
  _self.bc[0][0] = 0.0
  _self.bc[1][0] = 0.0


 def neumann_condition(_self, _neumann_pts):
  _self.bc_neumann = np.zeros([_self.npoints,1], dtype = float) 
  _self.neumann_pts = _neumann_pts 
 
  for i in range(0, len(_self.neumann_pts)):
   line = _self.neumann_pts[i][0] - 1
   v1 = _self.neumann_pts[i][1] - 1

   _self.bc_neumann[v1] += _self.bc[line]


 def dirichlet_condition(_self, _dirichlet_pts):
  _self.bc_dirichlet = np.zeros([_self.npoints,1], dtype = float) 
  _self.ibc = [] 
  _self.bc_1 = np.zeros([_self.npoints,1], dtype = float) #For scipy array solve
  _self.dirichlet_pts = _dirichlet_pts
 

  for i in range(0, len(_self.dirichlet_pts)):
   line = _self.dirichlet_pts[i][0] - 1
   v1 = _self.dirichlet_pts[i][1] - 1

   _self.bc_1[v1] = _self.bc[line]

   _self.bc_neumann[v1] = 0.0 #Dirichlet condition is preferential

   _self.ibc.append(v1)
   
  _self.ibc = np.unique(_self.ibc)


 def gaussian_elimination(_self, _LHS0, _neighbors_nodes):
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.bc_2 = np.ones([_self.npoints,1], dtype = float) 
  _self.neighbors_nodes = _neighbors_nodes

  for mm in _self.ibc:
   for nn in _self.neighbors_nodes[mm]:
    _self.bc_dirichlet[nn] -= float(_self.LHS[nn,mm]*_self.bc_1[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.bc_dirichlet[mm] = _self.bc_1[mm]
   _self.bc_2[mm] = 0.0
 

 def initial_condition(_self):
  _self.c = np.copy(_self.bc_dirichlet)
  _self.vx = np.zeros([_self.npoints,1], dtype = float)

  for i in range(0,_self.npoints):
   _self.vx[i] = 1.0


  for i in range(0, _self.npoints):
   x = _self.x[i]

   if x >= 2 and x <= 4:
    alpha = (x - 2)/(3 - 2)
    _self.c[i] = np.sin(alpha*(np.pi/2))






class NS2D:


 def __init__(_self, _numPhysical, _numNodes, _numVerts):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.numVerts = _numVerts


 def gaussianElimination(_self, _LHS0, _dirichletNodesVx, _dirichletNodesVy, _dirichletNodesPressure, _neighborsNodes, _neighborsNodesPressure, _aux1BCVx, _aux1BCVy, _aux1BCPressure):
  _self.dirichletNodesVx = _dirichletNodesVx 
  _self.dirichletNodesVy = _dirichletNodesVy 
  _self.dirichletNodesPressure = _dirichletNodesPressure
  _self.aux1BCVx = _aux1BCVx
  _self.aux1BCVy = _aux1BCVy
  _self.aux1BCPressure = _aux1BCPressure
  #_self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.LHS = np.copy(_LHS0)
  _self.neighborsNodes = _neighborsNodes
  _self.neighborsNodesPressure = _neighborsNodesPressure
  _self.dirichletVector = np.zeros([2*_self.numNodes + _self.numVerts,1], dtype = float) 
  _self.aux2BC = np.ones([2*_self.numNodes + _self.numVerts,1], dtype = float) 



  # Gaussian elimination for vx
  for mm in _self.dirichletNodesVx:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BCVx[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0

   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BCVx[mm]
   _self.aux2BC[mm] = 0.0




  # Gaussian elimination for vy
  for mm in _self.dirichletNodesVy:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn + _self.numNodes] -= float(_self.LHS[nn + _self.numNodes,mm + _self.numNodes]*_self.aux1BCVy[mm])
    _self.LHS[nn + _self.numNodes,mm + _self.numNodes] = 0.0
    _self.LHS[mm + _self.numNodes,nn + _self.numNodes] = 0.0
   
   _self.LHS[mm + _self.numNodes,mm + _self.numNodes] = 1.0
   _self.dirichletVector[mm + _self.numNodes] = _self.aux1BCVy[mm]
   _self.aux2BC[mm + _self.numNodes] = 0.0





  # Gaussian elimination for pressure
  for mm in _self.dirichletNodesPressure:
   for nn in _self.neighborsNodesPressure[mm]:
    _self.dirichletVector[nn + 2*_self.numNodes] -= float(_self.LHS[nn + 2*_self.numNodes,mm + 2*_self.numNodes]*_self.aux1BCPressure[mm])
    _self.LHS[nn + 2*_self.numNodes,mm + 2*_self.numNodes] = 0.0
    _self.LHS[mm + 2*_self.numNodes,nn + 2*_self.numNodes] = 0.0
   
   _self.LHS[mm + 2*_self.numNodes,mm + 2*_self.numNodes] = 1.0
   _self.dirichletVector[mm + 2*_self.numNodes] = _self.aux1BCPressure[mm]
   _self.aux2BC[mm + 2*_self.numNodes] = 0.0


class Cavity:

 # ------------------------------------------------------------------------------------------------------
 # Use:

 # # Applying vx condition
 # condition_xvelocity = bc_apply.Half_Cavity(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_xvelocity.neumann_condition(mesh.neumann_edges[1])
 # condition_xvelocity.dirichlet_condition(mesh.dirichlet_pts[1])
 # condition_xvelocity.gaussian_elimination(LHS_vx0,mesh.neighbors_nodes)
 # vorticity_ibc = condition_xvelocity.ibc

 # # Applying vy condition
 # condition_yvelocity = bc_apply.Half_Cavity(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_yvelocity.neumann_condition(mesh.neumann_edges[2])
 # condition_yvelocity.dirichlet_condition(mesh.dirichlet_pts[2])
 # condition_yvelocity.gaussian_elimination(LHS_vy0,mesh.neighbors_nodes)

 # # Applying psi condition
 # condition_streamfunction = bc_apply.Half_Cavity(mesh.nphysical,mesh.npoints,mesh.x,mesh.y)
 # condition_streamfunction.streamfunction_condition(mesh.dirichlet_pts[3],LHS_psi0,mesh.neighbors_nodes)
 # ------------------------------------------------------------------------------------------------------


 def __init__(_self, _numPhysical, _numNodes, _x, _y):
  _self.numPhysical = _numPhysical
  _self.numNodes = _numNodes
  _self.x = _x
  _self.y = _y
  _self.wallVelocity = 1.0
  _self.benchmark_problem = 'linear Cavity'


 def xVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 2 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Moving Wall
   elif line == 3:
    _self.aux1BC[v1] = _self.wallVelocity
    _self.aux1BC[v2] = _self.wallVelocity

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vx
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



 def yVelocityCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.lil_matrix.copy(_LHS0)
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Noslip 
   if line == 1 or line == 2 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Moving Wall
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for vy
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def streamFunctionCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

 # Dirichlet condition
  for i in range(0, len(_self.boundaryEdges)):
   line = _self.boundaryEdges[i][0]
   v1 = _self.boundaryEdges[i][1] - 1
   v2 = _self.boundaryEdges[i][2] - 1

   # Bottom Line
   # psi_bottom can be any value. Because, important is psi_top - psi_bottom.
   # In this case, psi_bottom is zero
   if line == 1:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0
 
    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)

   # Top Line
   # Ref: Batchelor 1967 pag. 76 eq. 2.2.8
   # psi_top is also zero, because the volume mass flux is null
   elif line == 3:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)


   # Right and Left lines
   # psi is also zero, because the volume mass flux is null
   elif line == 2 or line == 4:
    _self.aux1BC[v1] = 0.0
    _self.aux1BC[v2] = 0.0

    _self.dirichletNodes.append(v1)
    _self.dirichletNodes.append(v2)


  _self.dirichletNodes = np.unique(_self.dirichletNodes)


  # Gaussian elimination for psi
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 


 def pressureCondition(_self, _boundaryEdges, _LHS0, _neighborsNodes):
  _self.dirichletVector = np.zeros([_self.numNodes,1], dtype = float) 
  _self.dirichletNodes = [] 
  _self.aux1BC = np.zeros([_self.numNodes,1], dtype = float) #For scipy array solve
  _self.aux2BC = np.ones([_self.numNodes,1], dtype = float) 
  _self.LHS = sps.csr_matrix.copy(_LHS0) #used csr matrix because LHS = lil_matrix + lil_matrix
  _self.boundaryEdges = _boundaryEdges
  _self.neighborsNodes = _neighborsNodes

  # Dirichlet condition
  _self.aux1BC[0] = 0.0 # node (0,0) null pressure
  _self.dirichletNodes.append(0)
  _self.dirichletNodes = np.unique(_self.dirichletNodes)

  # Gaussian elimination for pressure
  for mm in _self.dirichletNodes:
   for nn in _self.neighborsNodes[mm]:
    _self.dirichletVector[nn] -= float(_self.LHS[nn,mm]*_self.aux1BC[mm])
    _self.LHS[nn,mm] = 0.0
    _self.LHS[mm,nn] = 0.0
   
   _self.LHS[mm,mm] = 1.0
   _self.dirichletVector[mm] = _self.aux1BC[mm]
   _self.aux2BC[mm] = 0.0
 



