****************************************************************************
units are: mm, ms, MPa, N, g

****************************************************************************
The difference between Wall and NoWall is that Wall has the wall penalty
contact BC of Tahoe, whereas NoWall is just directly the displacement BC.

We thought the wall penalty BC would allow calculation of nodal reaction
forces, but this appears not to be the case for this element.  See comment
below.

****************************************************************************
to run tahoe:

tahoe_s -f QuarterCyl-Drained-Displ-NoWall.xml >& outscreen-drained.txt &

wait until finished

view results in paraview
QuarterCyl-Drained-Displ-NoWall.io0.exo

To run again, delete outscreen-drained.txt file first.

to run q8p8/parallel, type

(make sure appropriate openmpi modules are loaded first)

mpirun -np 4 tahoe_p -f QuarterCyl-Drained-Displ.xml >& outscreen-drained.txt &

****************************************************************************
To extract results for plotting in Matlab, type

translate < outnodes-QuarterCyl-Drained-Displ-NoWall

which will generate a .txt file that can be imported for plotting purposes.
If you want to try your hand at running translate, then just type "translate"
and follow the prompts.  The outnodes* files contain the result of entering
responses to the prompts.

to determine the node number at the centerline bottom of the mesh, use
paraview 3.12, Selection Inspector under View, Create Selection, Global Node IDs, POINT, 
Invert Selection, Point Label, Visible, GlobalNodeID

Refer to the nodes_compare.m file for more information on what data are
plotted.

*************************************************************************
For the output results, "Kf-air" implies that the bulk modulus for the pore
"fluid" is for air

****************************************************************************
to extract axial force versus time curve from .exo file:

Note: it appears that the solid_fluid_mix poroelastic element in Tahoe does
not output nodal reaction forces at this time; needs to be added to source
code. But the following commands are provided anyway to show how it would be
done.

algebra QuarterCyl-Drained-Displ-NoWall.io1.exo R.exo

within algebra:
R=sum(F_D_Z)
save all
exit

blot R.exo
within blot:
tplot
typlot
R
pl
neutral
quit

this generates a R.xmgr text file that must be edited to remove extra lines for
plotting. After editing, rename the file

mv R.xmgr R.txt

then use gnuplot:

gnuplot plot_compare_single.gnu

edit plot_compare_single.gnu file to plot against Abaqus result.

Or import into Matlab for plotting.

****************************************************************************

