function [hconv]=hDB(h , p, speed, Dhy)
% hconv given by 0.023 RE^0.8 Pr^0.4 k / Dhy
% 
T    = tliq(h,p);
rhol = rholiq(h,p);
G    = rhol * speed;
mu   = muliq(T);
% 
% 
Re = G*Dhy/mu;
% 
Cp = cpliq(p,T);
kl = kliq(h);
%
Pr = mu * Cp / kl;
% 
hconv = 0.023 * Re^0.8 * Pr^0.4 * kl / Dhy;