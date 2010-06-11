function [vec] = solve_ss_hy( MSH, therm, hydro, vec, flxth, NLNewton, PreCond );

%%%%%%%%%%%%%%%%%%%%%%%%%%% WORK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note:  sum(flxth)*2*pi/sum(Pow)*therm.r_clad*MSH.dz*264/100 = 1.

if(MSH.dim==1)
    nx = 1;
else
    nx = MSH.nx;
end

Tin = hydro.T_in; % make radial dep later
hin = hydro.h_in;

if(NLNewton==true)
    H=spalloc(MSH.nz+1,MSH.nz+1,2*MSH.nz+1);
    H=speye(MSH.nz+1);
    for i=2:MSH.nz+1,
        H(i,i-1)=-1;
    end
    rhs=zeros(MSH.nz+1,1);
    rhs(1)=hin;
end

% radial loop
for ix=1:nx
    
    dec = (ix-1)*(MSH.nz+1) ;
    if(NLNewton==false)
        vec(dec+1) = hin;
        % axial loop
        for iz=1:MSH.nz
            flx = flxth(iz,ix);
            vec(dec+iz+1) = vec(dec+iz) + flx * ...
                2*pi*therm.r_clad*MSH.dz/100*therm.pindens*therm.ass_pitch*100 ...
                /hydro.mass_flow_rate;
        end
    else
        if(PreCond)
            vec(dec+1:dec+MSH.nz+1)=H\vec(dec+1:dec+MSH.nz+1);
        else
            for iz=1:MSH.nz
                flx = flxth(iz,ix);
                rhs(iz+1) =  flx * ...
                    2*pi*therm.r_clad*MSH.dz/100*therm.pindens*therm.ass_pitch*100 ...
                    /hydro.mass_flow_rate;
            end
            vec(dec+1:dec+MSH.nz+1)=H*vec(dec+1:dec+MSH.nz+1)-rhs;
        end
    end
end
