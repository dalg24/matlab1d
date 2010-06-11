function [vec]=solve_ss_hc( MSH, therm, hydro, vec, ent, Pow, Tmod, NLNewton, PreCond, x_kry);

if(PreCond & nargin==9)
    error('solve_ss_hc: with PreCond, x_hrylov must be given as 10th argument');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% NUMERICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           --------------------
% max. number of iterations to solve the non linear fuel conduction pb
therm.max_it_tq = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%% WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% therm.pindens: fuel pin density (npins/radial cm)
if(MSH.dim==1)
    coef = MSH.dz/100 * therm.pinsurf * therm.pindens * (therm.ass_pitch * 100) ;
    nx = 1;
% else
%     coef = MSH.dz/100 * therm.pinsurf * therm.pindens * MSH.dx ;
%     nx = MSH.nx;
end

dim_tq = MSH.radial.Ntot;
% radial loop (multi 1D channels)
for ix=1:nx
    dec     = (ix-1)*dim_tq*MSH.nz ;
    dec_ent = (ix-1)*(MSH.nz+1);
    % axial loop
    for iz=1:MSH.nz

        powdens = Pow(iz,ix) / coef ; % in W/m^3
        Tmod_ = Tmod(iz,ix);
        
        ib = dec + (iz-1)*dim_tq + 1;
        ie = dec + iz*dim_tq ;

        ib2 = iz   + dec_ent;
        ie2 = iz+1 + dec_ent;
        entm = 0.5*( ent(ib2)+ent(ie2) );
        hc = hDB( entm, hydro.pressure, hydro.speed, hydro.Dhy);

        if(NLNewton==false)
            Tf_old=vec(ib:ie);
            for itq=1:therm.max_it_tq,
                [ TQ, rhs_TQ ] = build_hc_ss( Tf_old, Tmod_, hc, ...
                    powdens, therm, MSH.radial);
                Tf = TQ\rhs_TQ;
                err = max( abs( Tf-Tf_old ) ) ;
                Tf_old = Tf;
                %             disp(sprintf('remaining error %g',err));
                if(err<1E-2)
                    %                 disp(sprintf('pb CVed in %i iterations',itq));
                    break
                elseif(itq==therm.max_it_tq)
                    disp(sprintf('remaining error %g',err));
                    disp(sprintf('%i iterations',itq));
                    disp('not enough iterations, pb not CVed in solve_ss_hc')
                end
            end
           vec(ib:ie) = Tf;
        else         
            [ TQ, rhs_TQ ] = build_hc_ss( vec(ib:ie), Tmod_, hc,...
                powdens, therm, MSH.radial);
            if(PreCond)
                vec(ib:ie) = TQ\x_kry(ib:ie) ;
            else
                vec(ib:ie) = TQ*vec(ib:ie) - rhs_TQ;
            end
        end
        
    end % axial loop
end % radial loop
