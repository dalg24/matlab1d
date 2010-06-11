function [vec] = init_vec( MSH, G, T_in, h_in, NPAR)

vec = ones(NPAR.siz,1);

skip_neu = NPAR.neu;

% set all Tfuel to T_in of the moderator + 100C
ib = skip_neu +1 ;
ie = skip_neu +  MSH.nz*MSH.radial.Ntot ;
vec(ib:ie) = T_in + 200;

% set all enth to h_in
ib = skip_neu + MSH.nz*MSH.radial.Ntot +1 ;
ie = NPAR.siz ;
vec(ib:ie) = h_in;

% % Teff = ones(MSH.nz,nx)*hydro.T_in+100;
% % Rmod = ones(MSH.nz,nx)*hydro.rho_in;
% % Tmod = ones(MSH.nz,nx)*hydro.T_in;
% % 
% % dim_tq = MSH.radial.Nfuel + MSH.radial.Nclad;
% % Tdist = zeros(dim_tq,MSH.nz,nx);
% % % Temp = ones(MSH.radial.Ntot,MSH.nz,MSH.nx)*hydro.T_in;

