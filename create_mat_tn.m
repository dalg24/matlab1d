function [A,ixs]=create_mat_tn(MSH,XS,NPAR,BC,Teff,Rmod,time,IVV);

% # of energy groups
n = NPAR.dof;
G = XS.G;
L = XS.nprec;

% memory allocation
siz = n*(G+L);
% Pp-M
if(MSH.dim==1)
    nnz = G^2*(n*3);
else
    nnz = G^2*(n*9);
end
% X, delayed n prod
nnz = nnz + n*L*G;
% B, delay precursor prod
if(MSH.dim==1)
    nnz = nnz + G*L*(n*3);
else
    nnz = nnz + G*L*(n*9);
end
% Lambda
nnz = nnz + n*L;

A=spalloc(siz,siz,nnz);
A(:,:)=0;

% interpolate XS
% loc=[18,1]*1;
loc=XS.loc;
[ixs]=InterpoleXS_Driver(MSH,NPAR,XS,Teff,Rmod,loc,time);

% assemble leakage
for g=1:G,
    ib=(g-1)*n+1;
    jb=ib;
    [A]=assemble_stiffness( A, ib, jb, -ixs.diff(:,g), MSH, NPAR );
end
% assemble removal
for g=1:G,
    ib=(g-1)*n+1;
    jb=ib;
    [A]=assemble_mass( A ,ib ,jb , -ixs.rem(:,g), MSH, NPAR);
end
if(G>1)
    % assemble downscattering
    ib=n+1;
    jb=1;
    [A]=assemble_mass( A, ib, jb, +ixs.scatt(:), MSH, NPAR);
end

% assemble prompt production
for ig=1:G,
    for jg=1:G,
        ib=(ig-1)*n+1;
        jb=(jg-1)*n+1;
        % only one set of beta is used !!!!
        data=XS.chi(1,ig)*(1-XS.btot(1)).*ixs.nuf(:,jg);  
        [A]=assemble_mass( A, ib, jb, data, MSH, NPAR);
    end
end

% assemble X part (delayed neutron production)
for i=1:L,
    for g=1:G,
        ib=(g-1)*n+1;
        jb=G*n+(i-1)*n+1;
        data=XS.chi(1,g)*XS.lambda(1,i)*ones(NPAR.nel,1);
        [A]=assemble_mass( A, ib, jb, data, MSH, NPAR, 1);
    end
end

% assemble B part (precursor production)
for g=1:G,
    for i=1:L,
        ib=G*n+(i-1)*n+1;
        jb=(g-1)*n+1;
        data=XS.b(1,i).*ixs.nuf(:,g);
        [A]=assemble_mass( A, ib, jb, data, MSH, NPAR);
    end
end

% assemble precursor removal
for i=1:L,
    ib=G*n+(i-1)*n+1;
    jb=ib;
    [A]=assemble_mass( A, ib, jb, -XS.lambda(1,i)*ones(NPAR.nel,1), MSH, NPAR, 1); 
end

if(nargin>7)
    A=IVV*A;
end

% % apply BC
if(MSH.dim==1)
    bc=BC.nodes(1,1:2);
else
    bc=[BC.nodes(:,1) ; BC.nodes(:,4) ];
%     bc=[BC.nodes(:,1) ; BC.nodes(:,2); BC.nodes(:,3); BC.nodes(:,4) ];
end
[n1 n2]=size(bc);
bcnodes=reshape(bc,n1*n2,1);
a=find(bc>0);

for i=1:length(a),               % loop fo the number of constraints
    id=bcnodes(a(i));         % extract the dof of a constraint
    for ig=1:G+L,
%         id=id+(ig-1)*NPAR.dof;
        id=bcnodes(a(i))+(ig-1)*NPAR.dof;
        A(id,:)=0.0;       % set all the id-th row to zero
        A(:,id)=0.0;       % set all the id-th column to zero
        A(id,id)=1;        % set the id-th diagonal to unity
    end
end
