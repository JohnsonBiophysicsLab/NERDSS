%fit to an exponential growth, after a lag time tau.
%in beta(2) returns the time-scale of growth, units of s. 1/beta(2) is growth rate k.
%offset in y is beta(4).
%extent of growth is beta(1)-beta(4). This is limited in simulations by the number of clathrin that can fit on the membrane. 
%Fit against y, copy numbers of clathrin on the membrane. 
function y=myexp_lag(x, beta1, beta2, beta3, beta4)

y=beta4+heaviside(x-beta3).*beta1.*(1-exp(-(x-beta3)/beta2));

end
