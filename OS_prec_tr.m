function [x_krylov]=OS3(MSH, NPAR, XS, therm, hydro, BC, x_krylov, vec, time, dt, IV);
% vec   = last newton
% x_krly= from gmres
%%%%%%%%% build neutronics system

[Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);
[Anp1,ixs]   = create_mat_tn(MSH, XS, NPAR, BC, Teff, Rmod,time,IV);

strid = NPAR.prec;
x_krylov(1:strid) = ( speye(strid)-dt/2*Anp1 ) \ x_krylov(1:strid);

% % % P contains the cell averaged powers per axial cell and radial channel
% for block jacobi precond, the Pow is useless b/c the rhs is the input x_krylov!
[Pow,Plevel] = comp_power( ixs ,MSH, NPAR, x_krylov(1:strid), XS.G );
[Pow,Plevel] = norma_power( Pow, Plevel*therm.Watt0, false );

% % %%%%%%%%% extract heat flux
[flxth]=extract_flxth( MSH, hydro, vec(strid+1:NPAR.ther), Tmod, vec(NPAR.ther+1:NPAR.siz) );

%%%%%%%%% heat conduction res
bidon_pown = 0;
bidon_xn   = 0;
[x_krylov(NPAR.prec+1:NPAR.ther)] = solve_tr_hc( MSH, therm, hydro,...
    vec(NPAR.prec+1:NPAR.ther), bidon_xn, dt, vec(NPAR.ther+1:NPAR.siz), ...
    bidon_xn, Pow, bidon_pown, Tmod, NPAR.MatFree, true, ...
    x_krylov(NPAR.prec+1:NPAR.ther), NPAR.method );

%%%%%%%% solve tr enthalpies
bidon_flx    = 0;
bidon_flxn   = 0;

[x_krylov(NPAR.ther+1:NPAR.siz)] = solve_tr_hy( MSH, therm, hydro, ...
    x_krylov(NPAR.ther+1:NPAR.siz), bidon_xn, dt, bidon_flx, bidon_flxn, ...
    NPAR.MatFree, true, NPAR.method );
