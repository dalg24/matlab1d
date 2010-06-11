function [NPAR,BC]=numpar(MSH, groups)

% nonlinear solves done using Newton and not Operator-Split
% the variable below can be modified
NPAR.NLNewton = false; 

NPAR.MatFree = true; 
NPAR.PreCond = false;

% check
if((NPAR.MatFree) && (NPAR.NLNewton==false))
    disp('numpar:: It makes no sense to have MatFree=true and NLNewton=false');
    NPAR.MatFree = false;
    NPAR.PreCond = false;
end

% maximum # of steady-state coupled N-TH iterations to solve the steady
% state problem
NPAR.max_ss_iter = 150;

NPAR.MaxNewton=50;

if(NPAR.NLNewton)
    NPAR.renorm=false;
else
%     NPAR.max_it_tq = 10;
    NPAR.renorm=true;
end

if(MSH.dim==1)

    NPAR.m = eye(2);
    NPAR.k = [1 -1; -1 1]/2;
    NPAR.f = [1;1];

    % number of elements
    NPAR.nel = MSH.nz;
    % number of degrees of freedom
    NPAR.dof = MSH.nz+1;
    % build connectivity array for FEM
    NPAR.gn = zeros(NPAR.nel,2);
    for iel=1:NPAR.nel
        NPAR.gn(iel,:)=[iel,iel+1];
    end
    BC.nodes(1,1)=1;
    BC.nodes(1,2)=NPAR.dof;

end

% total length of the vector of unknowns
% dof*ng  : Flux
% nz*Ntot : Tfuel
% nz+1    : Enthalpy
NPAR.siz = (groups*NPAR.dof) + (MSH.nz*MSH.radial.Ntot) + (MSH.nz+1);

% index of the last neutronic variable in the solution vector
NPAR.neu  = groups*NPAR.dof ;
% index of the last heat conduction variable in the solution vector
NPAR.ther = groups*NPAR.dof + MSH.nz*MSH.radial.Ntot ;

% scaling factors to help with matrix conditioning
NPAR.scaling = [1;1e3;1e4];
NPAR.scaling = ones(1,3);

% damping for Newton direction or OS
NPAR.damp=0.75;

%NPAR.log_epsi=false;

