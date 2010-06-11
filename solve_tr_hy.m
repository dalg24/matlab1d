function [vec] = solve_tr_hy( MSH, therm, hydro, vec, vecn, deltat, ...
    flxth, flxthn, MatFree, PreCond, method )
%%%%%%%%%%%%%%%%%%%%%%%%%%% WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note:  sum(flxth)*2*pi/sum(Pow)*therm.r_clad*MSH.dz*264/100 = 1.

if(MSH.dim==1)
    nx = 1;
else
    nx = MSH.nx;
end

Tin = hydro.T_in; % make radial dep later
hin = hydro.h_in;

H = spalloc(MSH.nz+1,MSH.nz+1,2*MSH.nz+1);
H = speye(MSH.nz+1);
for i=2:MSH.nz+1,
    H(i,i-1) = -1;
end
rhs=zeros(MSH.nz+1,1);
rhs(1)=hin * hydro.speed/(MSH.dz/100);
if( method == 'CN' )
    rhsn=rhs;
end

H=H*hydro.speed/(MSH.dz/100);

% max_iter_hy = 10;

% therm.ass_pitch*100: assy pitch in cm
% thermiq.pindens = Ncray/fa_size; % pin per cm of radial direction
% --> therm.pindens*therm.ass_pitch*100 : nbr of pins
% 2*pi*therm.r_clad*MSH.dz/100: lateral surface, in m^2
% radial loop
for ix=1:nx

    dec = (ix-1)*(MSH.nz+1) ;
    if(MatFree==false)
        %         for iter_hy=1:max_iter_hy,

        %         aux = vec(dec+1:dec+MSH.nz+1);

        % axial loop
        for iz=1:MSH.nz
            flx = flxth(iz,ix);
            rhs(iz+1) =  flx * ...
                2*pi*therm.r_clad*therm.pindens*therm.ass_pitch*100 ...
                /hydro.mass_flow_rate*hydro.speed;
            if( method == 'CN' )
                flx = flxthn(iz,ix);
                rhsn(iz+1) =  flx * ...
                    2*pi*therm.r_clad*therm.pindens*therm.ass_pitch*100 ...
                    /hydro.mass_flow_rate*hydro.speed;
            end
        end
        if( method =='IE' )
            vec(dec+1:dec+MSH.nz+1) = (speye(MSH.nz+1)+deltat*H) ...
                \ (  vecn(dec+1:dec+MSH.nz+1) +deltat*rhs ) ;
        elseif( method == 'CN' )
            vec(dec+1:dec+MSH.nz+1) = ( speye(MSH.nz+1)+deltat/2*H ) ...
                \ ( (speye(MSH.nz+1)-deltat/2*H )* vecn(dec+1:dec+MSH.nz+1) +deltat*(rhs+rhsn)/2 ) ;
        end

        %         err=norm(aux-vec(dec+1:dec+MSH.nz+1),inf);
        %         if( err<1e-10)
        %             break
        %         elseif(iter_hy==max_iter_hy)
        %             disp(sprintf('remaining error %g',err));
        %             disp(sprintf('%i iterations',iter_hy));
        %             disp('not enough iterations, pb not CVed in solve_tr_hy')
        %         end
        %     end
    else
        % axial loop
        if(PreCond)
            % Id -dt/2*df/dy is done is comp_res_diif_tr
            % that's why this line is slightly diff than the one above
            % for the precond case
            vec(dec+1:dec+MSH.nz+1) = ( speye(MSH.nz+1) +deltat/2*H ) \ vec(dec+1:dec+MSH.nz+1);
        else
            for iz=1:MSH.nz
                flx = flxth(iz,ix);
                rhs(iz+1) =  flx * ...
                    2*pi*therm.r_clad*therm.pindens*therm.ass_pitch*100 ...
                    /hydro.mass_flow_rate*hydro.speed;
            end
            vec(dec+1:dec+MSH.nz+1) = -( H*vec(dec+1:dec+MSH.nz+1) - rhs );
        end
    end
end

