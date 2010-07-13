function [res,err2] = comp_residual(MSH,XS,NPAR,BC,therm,hydro,vec);

res = zeros(NPAR.siz,1);

if(~NPAR.NLNewton)
    error('NPAR.NLNewton should be true, comp_residual.m')
end

%%%%%%%%% neutronics residual
% extract from the solution vector the fuel effective temperature and
% the moderator density and temperature
[Teff]      = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod,Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec );

%%%%%%%%% build neutronics system
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod );

strid = NPAR.neu;
res(1:strid) = (M-P) * vec(1:strid);

% P contains the cell averaged powers per axial cell and radial channel
[Pow,Plevel(1)] = comp_power( ixs ,MSH, NPAR, vec(1:strid), XS.G );
% useful to give close to the correct amplitude
Pow=Pow*therm.Watt0;
% [Pow,Plevel(1)] = norma_power( Pow, therm.Watt0, false );

%%%%%%%%% extract heat flux
[flxth]=extract_flxth( MSH, hydro, vec(NPAR.neu+1:NPAR.ther), Tmod, vec(NPAR.ther+1:NPAR.siz) );

%%%%%%%%% heat conduction res
[res(NPAR.neu+1:NPAR.ther)] = solve_ss_hc( MSH, therm, hydro,...
    vec(NPAR.neu+1:NPAR.ther), vec(NPAR.ther+1:NPAR.siz), Pow, Tmod, NPAR.NLNewton, false );
res(NPAR.neu+1:NPAR.ther) = res(NPAR.neu+1:NPAR.ther) / NPAR.scaling(2) ;

% %%%%%%%%% solve ss enthalpies
[res(NPAR.ther+1:NPAR.siz)] = solve_ss_hy( MSH, therm, hydro, ...
    vec(NPAR.ther+1:NPAR.siz), flxth, NPAR.NLNewton, false );
res(NPAR.ther+1:NPAR.siz) = res(NPAR.ther+1:NPAR.siz) / NPAR.scaling(3) ;



% compute the maximum of the ratio of residual(i)/solution_vector(i)
if(nargout>1)

    % we remove for the list of variables the ones for which solution_vector(i)=0
    % that is, the Dirchilet BC for the fluxes
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

    % compute ratio of residual/vector
    err2=0;
    for q=1:NPAR.siz
        ind=sum(find(id==q));
        if(ind==0)
            if(abs(vec(q))>1e-10)
                err2=max(err2, abs(res(q)/vec(q)) );
            end
        end
    end
end
