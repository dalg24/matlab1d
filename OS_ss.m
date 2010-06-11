function [vec,keff,Plevel] = OS_ss(MSH, NPAR, XS, therm, hydro, BC, vec);

% called here to get Teff, Tmod, Rmod, and  ixs
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod);

%%%%%%%%% solve ss neutronics                           
[keff,vec(1:NPAR.neu)] = solve_ss_neutro(MSH,NPAR,XS,BC,therm,hydro,vec);
disp(sprintf('ss eigenvalue is %g',keff));
%     [fl]=re_ordering(MSH,v,NPAR.gn,XS.G,NPAR.dof);

% Pow contains the cell averaged powers per axial cell and radial channel
[Pow,Plevel(1)] = comp_power( ixs, MSH, NPAR, vec(1:NPAR.neu), XS.G );
[Pow,Plevel(1)] = norma_power( Pow, therm.Watt0, true );

% % figure(11)
% % plot(Pow)
%%%%%%%%% solve ss heat conduction
[vec(NPAR.neu+1:NPAR.ther)] = solve_ss_hc( MSH, therm, hydro,...
    vec(NPAR.neu+1:NPAR.ther), vec(NPAR.ther+1:NPAR.siz), Pow, Tmod, NPAR.MatFree, NPAR.PreCond );
%%%%%%%% extract heat flux
[flxth]=extract_flxth( MSH, hydro, vec(NPAR.neu+1:NPAR.ther), Tmod, vec(NPAR.ther+1:NPAR.siz) );
%%%%%%%% solve ss enthalpies
[vec(NPAR.ther+1:NPAR.siz)] = solve_ss_hy( MSH, therm, hydro, ...
    vec(NPAR.ther+1:NPAR.siz), flxth, NPAR.MatFree, NPAR.PreCond );

