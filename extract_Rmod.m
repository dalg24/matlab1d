function [Rmod,Tmod]=extract_Rmod( MSH, strid, pressure, vec)

%%%% build Rmod
Rmod = zeros(MSH.nz,1);
Tmod = zeros(MSH.nz,1);

for iz=1:MSH.nz
    ib=iz   + strid;
    ie=iz+1 + strid;
    entm = 0.5*(vec(ib)+vec(ie));
    Rmod(iz,1) = rholiq(entm,pressure);
    Tmod(iz,1) = tliq(entm,pressure);
end
