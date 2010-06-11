function [resf,err2,therm] = comp_residual_tr(MSH,XS,NPAR,BC,therm,hydro,xnp1,time,dt,IV,src);

resf = zeros(NPAR.siz,1);

if(~NPAR.MatFree)
    error('NPAR.MatFree should be true, comp_residual.m')
end

%%%%%%%%% neutronics res
[Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xnp1 );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xnp1);
[Anp1, ixs]  = create_mat_tn(MSH,XS,NPAR,BC,Teff,Rmod,time,IV);


strid = NPAR.prec;
ns=length(diag(IV));
if(strid~=ns)
    error('strid~=ns')
end
resf(1:strid) = Anp1(1:ns,1:ns) * xnp1(1:strid);

% P contains the cell averaged powers per axial cell and radial channel
Plevel = therm.Plevel;
it     = length(Plevel);
[Pow,Plevel(it+1)] = comp_power( ixs, MSH, NPAR, xnp1(1:NPAR.neu), XS.G );
[Pow,Plevel(it+1)] = norma_power( Pow, Plevel(it+1)*therm.Watt0, false );

%%%%%%%%% heat conduction res
bidon_xkry = 0;
bidon_pown = 0;
bidon_xn   = 0;
[resf(NPAR.prec+1:NPAR.ther)] = solve_tr_hc( MSH, therm, hydro,...
    xnp1(NPAR.prec+1:NPAR.ther),  bidon_xn, dt, xnp1(NPAR.ther+1:NPAR.siz), ...
    bidon_xn, Pow, bidon_pown, Tmod, NPAR.MatFree, false, bidon_xkry, NPAR.method );

%%%%%%%%% extract heat flux
[flxth]=extract_flxth( MSH, hydro, xnp1(NPAR.prec+1:NPAR.ther), Tmod, ...
                        xnp1(NPAR.ther+1:NPAR.siz) );
                    
% %%%%%%%%% solve enthalpies
bidon_flxn = 0;
[resf(NPAR.ther+1:NPAR.siz)] = solve_tr_hy( MSH, therm, hydro, ...
    xnp1(NPAR.ther+1:NPAR.siz), bidon_xn, dt, flxth, bidon_flxn, ...
    NPAR.MatFree, false, NPAR.method );

% apply CN
if( nargin==10 )
    resf = xnp1 + dt / 2 * resf ;
elseif( nargin == 11 )
    resf = xnp1 - dt / 2 * resf -src;
else
    error('nargin error');
end

%
if( nargout>1 )

    if(MSH.dim==1)
        bc=BC.nodes(1,1:2);
    else
        bc=[BC.nodes(:,1) ; BC.nodes(:,4) ];
        %     bc=[BC.nodes(:,1) ; BC.nodes(:,2); BC.nodes(:,3); BC.nodes(:,4) ];
    end
    [n1 n2]=size(bc);
    bcnodes=reshape(bc,n1*n2,1);
    a=find(bc>0);
    id0=bcnodes(a(:));         % extract the dof of a constraint
    id=id0;
    for ig=1:XS.G,
        aux=id0+(ig-1)*NPAR.dof;
        id=[id;aux];
    end

    err2=0;
    for q=1:NPAR.siz
        ind=sum(find(id==q));
        if(ind==0)
            if(abs(xnp1(q))>1e-10)
                err2=max(err2, abs(resf(q)/xnp1(q)) );
            end
        end
    end
end

if( nargout >2 )
    therm.Plevel = Plevel ; 
end