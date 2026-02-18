%create rotation matrix for angle (in rads) around the z axis. 
function[M]=rotateEulerZ(angle)
format long
s=sin(angle);
c=cos(angle);

M=[c -s 0 ;
    s c 0;
    0 0 1]

    