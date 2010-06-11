function [ TQ, rhs_TQ ] = build_hc_tr( Tf, Tmod, hc, powdens, therm, geom, MatFree);

[ TQ, rhs_TQ ] = build_hc_ss( Tf, Tmod, hc, powdens, therm, geom );

Rfuel = therm.r_fuel;
Rclad = therm.r_clad;
Nfuel = geom.Nfuel;
Nclad = geom.Nclad;
d = Nfuel+Nclad;
e = Rclad-Rfuel;

M = speye(d);

% first ring
i = 1;
S = (Rfuel/(Nfuel-1))^2;
M(i,i) = RhoFuel(Tf(i)) * CpFuel(Tf(i)) * S/8;
% 2 to Nf-1 rings
for i=2:geom.Nfuel-1,
    S = (Rfuel/(Nfuel-1))^2 *(i-1);
    M(i,i) = RhoFuel(Tf(i)) * CpFuel(Tf(i)) * S;
end
% last ring
i = Nfuel;
S = (Rfuel/(Nfuel-1))^2 *(Nfuel-1-1/4)/2;
M(i,i) = RhoFuel(Tf(i)) * CpFuel(Tf(i)) * S;


% first ring
i = 1+Nfuel;
S = (Rfuel + e/(Nclad-1)/2)^2 - Rfuel^2;
M(i,i) = RhoClad(Tf(i)) * CpClad(Tf(i)) * S/2;
% 2 to Nf-1 rings
for i=2:Nclad-1,
    ii=i+Nfuel;
    S = (Rfuel + e/(Nclad-1)*(2*i-1)/2)^2 - (Rfuel + e/(Nclad-1)*(2*i-3)/2)^2 ;
    M(i,i) = RhoClad(Tf(ii)) * CpClad(Tf(ii)) *S/2 ;
end
% last ring
i = Nfuel+Nclad;
S = (Rfuel+Rclad)^2 - (Rfuel + e/(Nclad-1)*(2*Nclad-1)/2)^2;
M(i,i) = RhoClad(Tf(i)) * CpClad(Tf(i)) * S/2;


rhs_TQ = M \ rhs_TQ ;
TQ     = M \ TQ ;
