function [res] = comp_res_diff_tr(x_krylov, MSH, XS, NPAR, BC, therm, hydro, vec, time, dt, IV, resn, res0, epsi);


if(NPAR.PreCond==true)
    [x_krylov]=OS_prec_tr(MSH, NPAR, XS, therm, hydro, BC, x_krylov, vec, time, dt, IV);
end

% % if(NPAR.log_epsi)
% %     nor=norm(x_krylov,2);
% %     len=length(vec);
% %     s=1e-6*norm(vec,1)/len;
% %     if(nor>1e-12)
% %         epsi_=1e-6 + s/nor;
% %     else
% %         epsi_=1e-6;
% %     end
% % %     disp(sprintf(' computing epsi_: %g, %g, %g',s,nor,epsi_));
% %     
% %     if(nor>1e-12)
% %         epsi_=sqrt( (1+norm(vec,2))*1e-6 ) / nor;
% %     else
% % %         epsi_=sqrt( (1+norm(vec,2))*1e-6 ) / norm
% %         epsi_=1e-6;
% %     end
% %     disp(sprintf(' computing epsi_: %g, %g, %g',s,nor,epsi_));
% %     
% % else
% %     epsi_=epsi;
% % end
% % 

res = ( comp_residual_tr(MSH,XS,NPAR,BC,therm,hydro,vec+epsi*x_krylov,time,dt,IV,resn) - res0 ) / epsi ;

