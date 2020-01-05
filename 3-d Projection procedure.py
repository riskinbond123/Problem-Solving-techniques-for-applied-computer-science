# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:08:03 2019

@author: risha
"""

import numpy as np
import matplotlib.pyplot as plt

def projection_matrix(n,v):
    return (np.matmul(v,n.transpose()) - (np.dot(n.transpose(),v))*np.identity(4)).transpose()

tetrahedron = np.array([[0,0,0,1],[1,0,0,1],[0,1,0,1],[1,1,1,1]]).transpose()

"""
The question is presented in cartesian co-ordinates and it has been converted to their
respective homogeneous co-ordinates. Points at infinity have their last co-ordinate 0.
 
"""
# q1) Perspective projection onto the viewplane−x+ 3y+ 2z−4 = 0 from the viewpoint(2,−1,1).

viewplane_q1 = np.array([[-1,3,2,-4]]).transpose()
viewpoint_q1 = np.array([[2,-1,1,1]]).transpose()
print(projection_matrix(viewplane_q1,viewpoint_q1))
print("\n")


tetra_1 = np.matmul(tetrahedron.transpose(),projection_matrix(viewplane_q1,viewpoint_q1))
A_1 = (tetra_1[0]/tetra_1[0][3])[:3]
B_1 = (tetra_1[1]/tetra_1[1][3])[:3]
C_1 = (tetra_1[2]/tetra_1[2][3])[:3]
D_1 = (tetra_1[3]/tetra_1[3][3])[:3]

print(A_1)
print(B_1)
print(C_1)
print(D_1)
print("\n")
ax_1 = plt.axes(projection='3d')
ax_1.scatter3D([0,1,0,1], [0,0,1,1],[0,0,0,1] , cmap='orange')
ax_1.scatter3D([-2.67,-1.5,-0.33,1], [1.33,2.5,1.33,1],[-1.333,-1.333,-0.167,1] , cmap='Blues')


# q2)Perspective projection onto the viewplane 5x−3z+2 = 0 from the viewpoint (1,4,−1)

viewplane_q2 = np.array([[5,0,-3,2]]).transpose()
viewpoint_q2 = np.array([[1,4,-1,1]]).transpose()
print(projection_matrix(viewplane_q2,viewpoint_q2))
print("\n")


tetra_2 = np.matmul(tetrahedron.transpose(),projection_matrix(viewplane_q2,viewpoint_q2))
A_2 = (tetra_2[0]/tetra_2[0][3])[:3]
B_2 = (tetra_2[1]/tetra_2[1][3])[:3]
C_2 = (tetra_2[2]/tetra_2[2][3])[:3]
D_2 = (tetra_2[3]/tetra_2[3][3])[:3]

print(A_2)
print(B_2)
print(C_2)
print(D_2)
print("\n")
ax_2 = plt.axes(projection='3d')
ax_2.scatter3D([0,1,0,1], [0,0,1,1],[0,0,0,1] , cmap='orange')
ax_2.scatter3D([-0.25,1,-0.25,1], [-1,-9.333,0.25,-1],[0.25,2.333,0.25,2.3333] , cmap='Blues')


# q3) Parallel projection onto the viewplane 2y+ 3z+ 4 = 0 in the direction of the vector(1,−2,3).

viewplane_q3 = np.array([[0,2,3,4]]).transpose()
viewpoint_q3 = np.array([[1,-2,3,0]]).transpose()
print(projection_matrix(viewplane_q3,viewpoint_q3))
print("\n")


tetra_3 = np.matmul(tetrahedron.transpose(),projection_matrix(viewplane_q3,viewpoint_q3))
A_3 = (tetra_3[0]/tetra_3[0][3])[:3]
B_3 = (tetra_3[1]/tetra_3[1][3])[:3]
C_3 = (tetra_3[2]/tetra_3[2][3])[:3]
D_3 = (tetra_3[3]/tetra_3[3][3])[:3]

print(A_3)
print(B_3)
print(C_3)
print(D_3)
print("\n")
ax_3 = plt.axes(projection='3d')
ax_3.scatter3D([0,1,0,1], [0,0,1,1],[0,0,0,1] , cmap='orange')
ax_3.scatter3D([-0.8,0.2,-1.2,-0.8], [1.6,1.6,3.4,4.6],[-2.4,-2.4,-3.6,-4.4] , cmap='Blues')


# q4) Parallel projection onto the viewplane 7x−8y+ 5 = 0 in the direction of the vector(0,4,9)

viewplane_q4 = np.array([[7,-8,0,5]]).transpose()
viewpoint_q4 = np.array([[0,4,9,0]]).transpose()
print(projection_matrix(viewplane_q4,viewpoint_q4))
print("\n")


tetra_4 = np.matmul(tetrahedron.transpose(),projection_matrix(viewplane_q4,viewpoint_q4))
A_4 = (tetra_4[0]/tetra_4[0][3])[:3]
B_4 = (tetra_4[1]/tetra_4[1][3])[:3]
C_4 = (tetra_4[2]/tetra_4[2][3])[:3]
D_4 = (tetra_4[3]/tetra_4[3][3])[:3]

print(A_4)
print(B_4)
print(C_4)
print(D_4)
print("\n")
ax_4 = plt.axes(projection='3d')
ax_4.scatter3D([0,1,0,1], [0,0,1,1],[0,0,0,1] , cmap='orange')
ax_4.scatter3D([0,1,0,1], [0.625,1.5,0.625,1.5],[1.40625,3.375,-0.84375,2.125] , cmap='Blues')

