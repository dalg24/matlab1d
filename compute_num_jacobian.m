function J=compute_num_jacobian( unpert_res, MSH, XS, NPAR, BC, therm, hydro, vec)

JJ=zeros(NPAR.siz,NPAR.siz);
myeps=1e-6;

for i=1:NPAR.siz
    % save previous value
    sav = vec(i);
    % perturb i-th component
    vec(i)=vec(i) + myeps;
    % evaluate perturbed residual
    [res] = comp_residual( MSH, XS, NPAR, BC, therm, hydro, vec );
    % store i-th column of jacobian matrix
    JJ(:,i)=(res-unpert_res)/myeps;
    % undo the perturbation
    vec(i) = sav;
end

% output in sparse format to save memory
J=sparse(JJ);

    

