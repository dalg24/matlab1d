function [A]=assemble_mass(A,ib,jb,data,MSH,NPAR,opt)

n_el=NPAR.nel;
% z2med = MSH.z2mat;

if(nargin>6)
    if(opt==1)
        if(MSH.dim==1)
            jac = 1/2;
            my_m = eye(2);
        else
            jac = 1/4;
            my_m = eye(4);
        end
    else
        error('cannot have 7 inputs and opt/=1')
    end
else
    if(MSH.dim==1)
        jac = MSH.dz / 2;
    else
        jac = MSH.dz * MSH.dx / 4;
    end
    my_m=NPAR.m;
end

for iel=1:n_el,
    gn  = NPAR.gn(iel,:);
    ign=ib-1+gn;
    jgn=jb-1+gn;
%     i_med=z2med(iel);
    A(ign(:),jgn(:)) = A(ign(:),jgn(:)) + jac*data(iel)*my_m;
end
