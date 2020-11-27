# O(n) Collision Detection

This code was developed as an assignment for [Animation and CGI Motion](https://www.edx.org/course/animation-and-cgi-motion-2), a graduate level course from Columbia University. 
This code begins by breaking the animation region into a fixed grid, and determining particles that are either in or on the edge of the grid. Within each grid, a sweep and prune method is performed to determine what particles, edges, and halfplanes could be colliding. The sweep and prune method first determines potential intersections along the x axis, then potential intersections along the y axis. This is achieved by determining the min and max extent of each object along the axis, and determining possible overlaps. If an object is overlapping in both the X and Y axis, then a collision is found. This method is O(n) because the particles are assumed to be in the same order during each timestep. There assumed positions along the X and Y axis are updated through a sorting algorithm to determined their new positions, then collisions are checked for. 

Future work for this would be using a dynamic grid, and accounting for continuous time collision detection. 
