function [res] = comp_res_diff( x_krylov, MSH, XS, NPAR, BC, therm, hydro, vec, src, res0, epsi )

if(NPAR.PreCond)
    [x_krylov,bidon]=OS_prec_ss(MSH, NPAR, XS, therm, hydro, BC, x_krylov, vec);
end

res = ( comp_residual( MSH,XS,NPAR,BC,therm,hydro,vec+epsi*x_krylov,src) - res0 ) / epsi ;

