function [A]=assemble_stiffness(A,ib,jb,data,MSH,NPAR)


n_el=NPAR.nel;
z2med = MSH.z2mat;

if(MSH.dim==1)
    jac = MSH.dz / 2;
    for iel=1:n_el,
        gn  = NPAR.gn(iel,:);
        ign = ib-1+gn;
        jgn = jb-1+gn;
%         i_med=z2med(iel);
        A(ign(:),jgn(:)) = A(ign(:),jgn(:)) +1/jac*data(iel)*NPAR.k ;
    end
else
    jac = MSH.dz * MSH.dx / 4;
    ratio = MSH.dz / MSH.dx;
    for iel=1:n_el,
        gn  = NPAR.gn(iel,:);
        ign = ib-1+gn;
        jgn = jb-1+gn;
%         i_med=z2med(iel);
        A(ign(:),jgn(:)) = A(ign(:),jgn(:)) + data(iel) * ...
                           ( ratio*NPAR.kx  + 1/ratio*NPAR.ky ) ;
    end
end

% for iel=1:n_el,
%
%     gn  = NPAR.gn(iel,:);
%
% %     iz  = mod(iel,MSH.nz);
% %     ix  = floor((iel-1)/nz)+1;
%
%
%     ign = ib-1+gn;
%     jgn = jb-1+gn;
%     jacob=dx/2;
%     i_med=z2med(iel);
%     A(ign(:),jgn(:)) = A(ign(:),jgn(:)) +1/jacob*data(i_med)*k ;
% end

% + jacob*data(i_med)*m