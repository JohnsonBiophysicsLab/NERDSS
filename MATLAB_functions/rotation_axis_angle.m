function[M]=rotation_axis_angle(theta, u)
format long
ux=u(1);
uy=u(2);
uz=u(3);
ct=cos(theta);
st=sin(theta);
m=1-ct;

M(1,1)=ct+ux^2*m;
M(1,2)=ux*uy*m-uz*st;
M(1,3)=ux*uz*m+uy*st;
M(2,1)=uy*ux*m+uz*st;
M(2,2)=ct+uy^2*m;
M(2,3)=uy*uz*m-ux*st;
M(3,1)=uz*ux*m-uy*st;
M(3,2)=uz*uy*m+ux*st;
M(3,3)=ct+uz^2*m;

