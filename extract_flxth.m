function [flxth]=extract_flxth( MSH, hydro, Tf, Tmod, ent )

%%%% build heat flux
flxth = zeros(MSH.nz,1);

% hc = hDB( hydro.h_in, hydro.pressure, hydro.speed, hydro.Dhy);
% flxth(iz,ix) = hc*(Tf(dim_tq)-Tmod_);

dim_tq = MSH.radial.Ntot ;

for iz=1:MSH.nz

    % compute the average enthalpy in axial cell 
    ib=iz   ;
    ie=iz+1 ;
    entm = 0.5*( ent(ib)+ent(ie) );
    % compute the heat exchange coefficient
    hc = hDB( entm, hydro.pressure, hydro.speed, hydro.Dhy);

    % retrieve the fuel element surface temperature
    % and use it to compute the heat flux
    ii = iz*dim_tq ;
    flxth(iz,1) = hc*( Tf(ii) - Tmod(iz,1) );

end
