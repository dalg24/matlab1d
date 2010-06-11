function [P,M,ixs]=create_mat_ss(MSH,XS,NPAR,BC,Teff,Rmod,loc,time);


if(nargin>6)
    [ixs]=InterpoleXS_Driver(MSH,NPAR,XS,Teff,Rmod,loc,time);
else
    [ixs]=InterpoleXS_Driver(MSH,NPAR,XS,Teff,Rmod);
end    

% # of energy groups
G=XS.G;

% memory allocation
nnz=G^2*(NPAR.dof*3);
P=spalloc(G*NPAR.dof,G*NPAR.dof,nnz);
P(:,:)=0;
M=spalloc(G*NPAR.dof,G*NPAR.dof,nnz);
M(:,:)=0;

% assemble leakage
n=NPAR.dof;
for g=1:G,
    ib = (g-1)*n+1;
    jb = ib;
    [M]=assemble_stiffness( M, ib, jb, ixs.diff(:,g), MSH, NPAR );
end
% assemble removal
for g=1:G,
    ib = (g-1)*n+1;
    jb = ib;
    [M]=assemble_mass( M ,ib ,jb , ixs.rem(:,g), MSH, NPAR);
end
if(G>1)
    % assemble downscattering
    ib = n+1;
    jb = 1;
    [M]=assemble_mass( M, ib, jb, -ixs.scatt(:), MSH, NPAR);
end

% assemble production
for ig=1:G,
    for jg=1:G,
        ib = (ig-1)*n+1;
        jb = (jg-1)*n+1;
        data=XS.chi(1,ig).*ixs.nuf(:,jg);
        [P]=assemble_mass( P, ib, jb, data, MSH, NPAR);
    end
end

bc=BC.nodes(1,1:2);
[n1 n2]=size(bc);
bcnodes=reshape(bc,n1*n2,1);
a=find(bc>0);

for i=1:length(a),               % loop fo the number of constraints
    id=bcnodes(a(i));         % extract the dof of a constraint
    for ig=1:G,
        id=id+(ig-1)*NPAR.dof;
        M(id,:)=0.0;       % set all the id-th row to zero
        M(:,id)=0.0;       % set all the id-th column to zero
        M(id,id)=1;        % set the id-th diagonal to unity
        P(id,:)=0.0;       % set all the id-th row to zero
    end
end

