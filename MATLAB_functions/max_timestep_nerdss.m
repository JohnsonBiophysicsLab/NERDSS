%Da in um2/s
%Db in um2/s
%rho in Ncopies/um3
%sigma in um
%OUTPUTS:
%deltat in s.
function[deltat3d]=max_timestep_nerdss(Da, Db, rhoa, rhob, sigma)


deltat3d=1/54/(Da+Db)* ((3/(4*pi*max(rhoa, rhob))+sigma^3)^(1/3)-sigma)^2;
deltat2d=1/36/(Da+Db)* ((1/(pi*max(rhoa, rhob))+sigma^2)^(1/2)-sigma)^2;
