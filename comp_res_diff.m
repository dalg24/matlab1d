function [res] = comp_res_diff( x_krylov, MSH, XS, NPAR, BC, therm, hydro, vec, res0, epsi )

if(NPAR.PreCond)
    [x_krylov]=OS_prec_ss(MSH, NPAR, XS, therm, hydro, BC, x_krylov, vec);
end

res = ( comp_residual( MSH,XS,NPAR,BC,therm,hydro,vec+epsi*x_krylov) - res0 ) / epsi ;

