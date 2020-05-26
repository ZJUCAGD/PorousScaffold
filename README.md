# PorousScaffold
Generating Porous scaffold based on B-spline solid using TPMS.
The projects are run on Windows operating system, and thes specified compiling platform is Visual Studio 2013.

# MarchingTetrahedrons
Openmesh, OpenGL and Eigen libraries are called.

`extractisosurface:`\
Generate the TPMS structures in the parametric domain by Marching Tetrahedral algorithm.\
`curvature2threshold:` \
Map the mean curvature of the solid model's surface to the threshold filed according to the TPMS, and a TDF file is generated.\
`porosityPereservation:` \
Preserve the porosity of the porous scaffold before and after mapping by the optimization method.

# ScaffoldGeneration
Qt5.0, OpenGl and Eigen libraries are called.

Read TDF files to generate porous scaffolds.
