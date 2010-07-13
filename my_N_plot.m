function my_N_plot(fid, clf_, hold_, MSH, NPAR, XS, BC, therm, hydro, vec )

% extract from the solution vector the fuel effective temperature and
% the moderator density and temperature
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);

figure(fid); 
% clear figure if requested
if(clf_), clf; end
subplot(3,1,1); plot(vec(1:NPAR.neu/2),'r-+'); 
% hold on if requested
if(hold_), hold on; end
subplot(3,1,2); plot(vec(NPAR.neu/2+1:NPAR.neu),'r-+'); 
% hold on if requested
if(hold_), hold on; end

[Teff]      = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod,Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec );
%%%%%%%%% build neutronics system
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod );
% P contains the cell averaged powers per axial cell and radial channel
[Pow,Plevel(1)] = comp_power( ixs ,MSH, NPAR, vec(1:NPAR.neu), XS.G );
% Pow=Pow*therm.Watt0;

subplot(3,1,3); plot(Pow,'r-+'); 

drawnow;

