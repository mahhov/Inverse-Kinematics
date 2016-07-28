# Inverse-Kinematics  
  
https://www.youtube.com/watch?v=HTdj-V_xGL8   
  
~May 2014  
    
controls:  
   
1) s -> toggles flat/smooth shading  
2) w -> toggles wireframe/solid  
3) q -> exits  
4) space -> begins/pauses/resumes (the program begins paused, press space initially when ready to begin)  
5) r -> restarts to initial configuration and resets path (if input file has been modified, will apply the changes)  
6) t -> toggles auto-rotates camera (horizontally)  
7) o -> randomizes arm position  
8) i -> enters "unending-mode" where random points will continuosly be added to the end of the path  
9) c -> toggles displaying/hiding the path  
10) arrow keys -> move path (horizontally)  
11), and . -> move path (vertically)  
12)shift + left mouse -> move path (horizontally)  
13)shift + right mouse -> move path (vertically)  
14)left mouse -> rotate camera  
15)middle mouse -> zoom camera  
16)right mouse -> pan camera  
17)+ (or =) and - -> speed up/down  
  
input file:  
  
1) first line -> first number determines how many arm pieces, following numbers determine length of each  
2) second line -> first number determines how "key" path points you have, second number determines how many path points to add in between consequtive "key" points (via linear interpolation)  
3) following lines -> x y z coordinates of each "key" path point  
  
algorithm used:  
  
1) compute Jacobian  
2) compute pesuedo inverse via J+ = J^T * (J * J^T)^-1  
3) find ratio of parameter changes via J+ * (goal - current)  
4) normalize and apply changes (if it improves position, otherwise use smaller step size)  
   
to run:  
  
1) either open the .sln in visual studio and run  
2) or on linux, type "make", then "./as1 input.txt .1 .0005 .015" (just an example)  
3) first argument is input file describing arms and path  
4) second argument is "tolerance" (how close our endpoint has to be to path), good value is around .1  
5) third argument is "step size" (rouch measure how much parameters change per jacobian), good value is around .0005.  
6) fourth argument is "default speed" (roughly, how fast animation occurs), good value is between .015 to .09. speed can be changed later by pressing +/-/=.  
