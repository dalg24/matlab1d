function J=compute_num_jacobian_eig( unpert_res, MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old)

JJ=zeros(NPAR.siz,NPAR.siz);
myeps=1e-6;

for i=1:NPAR.siz
    % save previous value
    sav = vec(i);
    % perturb i-th component
    vec(i)=vec(i) + myeps;
    % evaluate perturbed residual
    [res] = comp_residual_eig( MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old );
    % store i-th column of jacobian matrix
    JJ(:,i)=(res-unpert_res)/myeps;
    % undo the perturbation
    vec(i) = sav;
end

% output in sparse format to save memory
J=sparse(JJ);

    

