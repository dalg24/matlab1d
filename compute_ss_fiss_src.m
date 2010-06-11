function [fis_src,P,M,ixs]=compute_ss_fiss_src( MSH, NPAR, BC, XS, therm, hydro, vec)

% P: neutron production matrix
% M: neutron loss matrix
% ixs: inteprolated XS
% fis_src: fission fission

fis_src = zeros(NPAR.siz,1);

% extract from the solution vector the fuel effective temperature and
% the moderator density and temperature
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);

%%%%%%%%% build neutronics system
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod);

% compute the steady-state fission source
fis_src(1:NPAR.neu) = P*vec(1:NPAR.neu);
