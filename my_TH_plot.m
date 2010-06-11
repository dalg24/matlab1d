function my_TH_plot(fid, clf_, hold_, MSH, NPAR, therm, hydro, vec )

% extract from the solution vector the fuel effective temperature and
% the moderator density and temperature
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);

figure(fid); 
% clear figure if requested
if(clf_), clf; end
% plot Teff
subplot(2,1,1); plot(Teff,'r-+'); 
% hold on if requested
if(hold_), hold on; end
subplot(2,1,2); plot(Tmod,'r-+'); 
% hold on if requested
if(hold_), hold on; end

drawnow;

