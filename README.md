# 	Finite element method simulation of a simplified memristor model

## General description:
#### Author: Patrik Reizinger
#### Purpose: 
Created as the part of the course *Field Simulation with Finite Element Methods (BMEVIHVJV35)* at the __*Budapest University of Technology and Economics*__ __(BUTE)__.
For legal matters see the included *LICENSE* file, with the following __restrictions__: the use of this material for the above mentioned course or its successors is prohibited, except to study the project while seeking technical help, e.g how to use specific Matlab functions. But regarding the additional *content*, please do not use this material which was made available to help other enthusiasts gaining a better control above the *PDE Toolbox* of Matlab.


## memristor_pde.m
The __*memristor_pde.m*__ function creates a PDE model then prcesses it including the following steps:
- Geometry
- Mesh
- Boundary conditions
- PDE coefficient specification
- Solution

### Input arguments:
- **semi_a**: semiaxis of the ellipse (along the X-axis)
- **semi_b**: semiaxis of the ellipse (along the Y-axis)
- **boundary_offset**: offset of the material boundary (0 means the sigmoid is "centered" at 0)
- **plot_flag**: binary variable to switch the plot function on/off

### Output arguments:
- **results**: solution of the PDE
- **model**: PDE model

### Example use:
``` matlab
help memristor_pde % if you need help, it prints out both the functionality and the description of the parameters
[results, model] = memristor_pde(semi_a, semi_b, boundary_offset, plot_flag)
```

## process_results.m
The __*process_results.m*__ script makes the postporcessing and visualization on the solution(s) obtained with __*memristor_pde.m*__.

### Parameters
The following parameters can be specified to customize the created plots/animations:
- **plot_mode**:
    - __*0*__ - u-plot
    - __*1*__ - E + u-plot without c
    - __*2*__ - E + u-plot with c (D)
- **change_mode**:
    - __*0*__ - ellipse size will be changed
    - __*1*__ - state boundary will be changed
    - __*2*__ - ellipse and state boundary will be changed
- **en_3D**: while plotting, creates a 3D figure is set to nonzero
- **write_video**: flag to specify whether to write videos into file (.avi), target directory is ./videos
- **show_animation**: flag to specify whether to show the animations

### Example use:
``` matlab
help process_results % if you need help, it prints out both the functionality and the description of the parameters
process_results
```

