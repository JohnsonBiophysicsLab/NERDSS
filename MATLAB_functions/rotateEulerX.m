%create rotation matrix for angle (in rads) around the X axis. 
function[M]=rotateEulerX(angle)
format long
s=sin(angle);
c=cos(angle);

M=[1 0 0 ;
    0 c -s;
    0 s c]

    