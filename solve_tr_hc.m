function [vec]=solve_tr_hc( MSH, therm, hydro, vec, vecn, dt, ent, entn, Pow, ...
    Pown, Tmod, MatFree, PreCond, x_kry, method);

%%%%%%%%%%%%%%%%%%%%%%%%%%% NUMERICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           --------------------
% max. number of iterations to solve the non linear fuel conduction pb
therm.max_it_tq = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%% WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% therm.pindens: fuel pin density (npins/radial cm)
if(MSH.dim==1)
    coef = MSH.dz/100 * therm.pinsurf * therm.pindens * (therm.ass_pitch * 100) ;
    nx = 1;
else
    coef = MSH.dz/100 * therm.pinsurf * therm.pindens * MSH.dx ;
    nx = MSH.nx;
end
% if(MatFree==false)
%     flxth= zeros(MSH.nz,nx);
% end

dim_tq = MSH.radial.Ntot;
% radial loop
for ix=1:nx
    dec     = (ix-1)*dim_tq*MSH.nz ;
    dec_ent = (ix-1)*(MSH.nz+1);
    % axial loop
    for iz=1:MSH.nz

        powdens = Pow(iz,ix) / coef ; % in W/m^3
        if(MatFree==false)
            powdensn= Pown(iz,ix) / coef ; % in W/m^3
        end
        Tmod_ = Tmod(iz,ix);

        ib = dec + (iz-1)*dim_tq + 1;
        ie = dec + iz*dim_tq ;

        ib2 = iz   + dec_ent;
        ie2 = iz+1 + dec_ent;
        entm = 0.5*( ent(ib2)+ent(ie2) );
        if(MatFree==false)
            entmn = 0.5*( entn(ib2)+entn(ie2) );
            hcn   = hDB( entmn, hydro.pressure, hydro.speed, hydro.Dhy);
        end
        hc = hDB( entm, hydro.pressure, hydro.speed, hydro.Dhy);


        if(MatFree==false)

            [ TQn, rhs_TQn ] = build_hc_tr( vecn(ib:ie), Tmod_, hcn, ...
                powdensn, therm, MSH.radial, MatFree);
            Tf_old=vecn(ib:ie);
            for itq=1:therm.max_it_tq,
                [ TQ, rhs_TQ ] = build_hc_tr( Tf_old, Tmod_, hc, ...
                    powdens, therm, MSH.radial, MatFree);
                if(method=='IE')
                    rhs = dt*rhs_TQ + vecn(ib:ie);
                    Tf = ( speye(dim_tq) + dt*TQ)\rhs;
                elseif(method=='CN')
                    rhs = dt*(rhs_TQ+rhs_TQn)/2 + ( speye(dim_tq) -dt/2*TQn )*vecn(ib:ie);
                    Tf = ( speye(dim_tq) + dt/2*TQ)\rhs;
                end
                err = max( abs( Tf-Tf_old ) ) ;
                Tf_old = Tf;
                %             disp(sprintf('remaining error %g',err));
                if(err<1E-12)
                    %                 disp(sprintf('pb CVed in %i iterations',itq));
                    break
                elseif(itq==therm.max_it_tq)
                    disp(sprintf('remaining error %g',err));
                    disp(sprintf('%i iterations',itq));
                    disp('not enough iterations, pb not CVed in solve_tr_hc')
                end
            end

            vec(ib:ie) = Tf;

        else
            [ TQ, rhs_TQ ] = build_hc_tr( vec(ib:ie), Tmod_, hc,...
                powdens, therm, MSH.radial, MatFree);
            if(PreCond)
                vec(ib:ie) = ( speye(dim_tq) +dt/2*TQ ) \x_kry(ib:ie) ;
            else
                % Id -dt/2*df/dy is done is comp_res_diif_tr
                % that's why this line is slightly diff than the one above
                % for the precond case
                vec(ib:ie) = -( TQ*vec(ib:ie) - rhs_TQ );
            end
        end

    end % axial loop
end % radial loop