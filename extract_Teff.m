function [Teff]=extract_Teff(MSH,strid,rw,vec)

%%%% build Teff
Teff = zeros(MSH.nz,1);

for iz=1:MSH.nz
    ib=(iz-1)*MSH.radial.Ntot +1 + strid;
    ie=(iz  )*MSH.radial.Ntot    + strid;
    Tf = vec(ib:ie);
    Teff(iz) = rw(1)*Tf(1) + rw(2)*Tf(MSH.radial.Nfuel);
end
