function [A,ixs]=create_mat_invv(XS,MSH,NPAR,Teff,Rmod);

% interpolate XS
[ixs]=InterpoleXS_Driver(MSH,NPAR,XS,Teff,Rmod);

% # of energy groups
G=XS.G;
L=XS.nprec;

% memory allocation
n=NPAR.dof;
siz=n*(G+L);
if(MSH.dim==1)
    nnz=G^2*(n*3);
else
    nnz=G^2*(n*9);
end
A=spalloc(siz,siz,nnz);
A(:,:)=0;

% assemble 

for g=1:G,
    ib = (g-1)*n+1;
    jb = ib;
    [A]=assemble_mass( A, ib, jb, ixs.invvelo(:,g), MSH, NPAR );
end
for i=G*n+1:siz,
    A(i,i)=1;
end

% already compute the inverse of inv velocity
A=inv(A);