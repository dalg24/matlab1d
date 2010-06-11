function [p,P0]=comp_power(ixs,MSH,NPAR,v,G)

p  = zeros(NPAR.nel,1);

n_el=NPAR.nel;
z2med = MSH.z2mat;
for g=1:G,
    ib = (g-1)*NPAR.dof+1;
        jac = MSH.dz / 2;
    for iel=1:n_el,
        gn  = NPAR.gn(iel,:);
        ign=ib-1+gn;
%         i_med=z2med(iel);
        p(iel,1) = p(iel,1) + jac * ixs.fiss(iel,g) * ...
            dot( NPAR.f ,v(ign(:),1) );
    end
end

P0=sum(sum(p));

