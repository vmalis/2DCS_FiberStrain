2DCS_FiberStrain


Analysis pipeline

1.	RECON:: CS reconstruction;  input: raw k-space data, output: temporal frames of magnitude and phase images 
2.	Fiber tracking of medial Gastrocs fibers:: input: water suppressed images, output: start and end coordinates of fibers (3 fibers each in distal, middle and proximal regions).  These regions are approximately the top third, middle and distant third of the MG muscle.  See image and sample fibers in the medial gastrocnemius.
2a. Get the slope of the line connecting the distal, middle and proximal points. We will for the abstract just get the angle to the y-axis.
2b. Track the start and end points so that the fiber trajectory and fiber angle is known as a function of the cycle.
2c. How do we measure the aponeurosis separation? We can manually trace the aponeurosis around the middle 1/3 rd and fit a line to both and determine the separation.  Or in a simpler way, we could do this: Mark two points for the deep and superficial aponeurosis (mark the points where they are approximately parallel to each other and just a bit inside the muscle). Then treat these two lines like the manual fibers.  Only addition is to also find the distance between these two lines.  We can get 
3.	2D Strain rate (SR) and 2D strain (L): Process the velocity data to generate the 2D SR and strain data.
4.	ROIs are placed on the locations shown in the figure; these ROIs are tracked through the temporal frame. Get SR, L (fiber, in-plane, out-plane, angle of the fiber to the y-axis and angle between principle negative strain and the fiber). 
