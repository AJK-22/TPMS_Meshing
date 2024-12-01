Isosurface_2D.png: Initial 2D demo used to develop tools for moving along a particular isosurface and traversing
  in the perpendicular direction to generate an inflation layer at the surface.

Gmsh_export_2D.png: The same mesh as shown in the previous image, exported from MATLAB to Gmsh.

Gyroid_Surface_Flowlines.png: Progress made with the gyroid TPMS geometry to traverse along surface curvature lines
  while placing nodes at regular intervals, the first step towards obtaining a surface mesh.

Schwarz-P_Surface_Flowlines.png: Similar methodology applied to the Schwarz-P TPMS geometry, flowlines originating 
  from opposite edges are merged together to prevent an excessive number of nodes.

Schwarz-P_Surface_Wireframe.png: Initial extraction of a wireframe mesh from the surface nodes. Note the extreme reduction
  in triangular surface elements. 3D mesh in inflation region may then be obtained by moving along function gradient towards 
  another isolevel at each surface node. Current state of progress is extracting a completely watertight mesh and converting 
  it to the .msh format for evaluation of mesh statistics against exsiting meshing tools.
