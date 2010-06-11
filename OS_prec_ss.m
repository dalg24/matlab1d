function [x_krylov]=OS_prec_ss(MSH, NPAR, XS, therm, hydro, BC, x_krylov, vec);

% extract from the solution vector the fuel effective temperature and
% the moderator density and temperature
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);

%%%%%%%%% build neutronics system
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod);

strid = NPAR.neu;
x_krylov(1:strid) = M \ x_krylov(1:strid);

% % % P contains the cell averaged powers per axial cell and radial channel
% for block jacobi precond, the Pow is useless b/c the rhs is the input x_krylov!
[Pow,Plevel(1)] = comp_power( ixs ,MSH, NPAR, x_krylov(1:strid), XS.G );
[Pow,Plevel(1)] = norma_power( Pow, therm.Watt0, false );

%%%%%%%%% solve ss heat conduction
[x_krylov(NPAR.neu+1:NPAR.ther)] = solve_ss_hc( MSH, therm, hydro,...
    vec(NPAR.neu+1:NPAR.ther), vec(NPAR.ther+1:NPAR.siz), Pow, Tmod, ...
    true, NPAR.PreCond, x_krylov(NPAR.neu+1:NPAR.ther) );

% % %%%%%%%%% extract heat flux
% % [flxth]=extract_flxth( MSH, hydro, vec(NPAR.neu+1:NPAR.ther), Tmod, vec(NPAR.ther+1:NPAR.siz) );
% %  I believe this is useless for matfree precond
bidon_flx    = 0;

%%%%%%%% solve ss enthalpies
[x_krylov(NPAR.ther+1:NPAR.siz)] = solve_ss_hy( MSH, therm, hydro, ...
    x_krylov(NPAR.ther+1:NPAR.siz), bidon_flx, true, NPAR.PreCond );

