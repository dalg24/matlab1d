function [ixs]=InterpoleXS_Driver(MSH,NPAR,XS,Teff,Rmod,loc,time)

% interpolate
if(MSH.dim==1)
    nx = 1;
else
    nx = MSH.nx;
end
e2ij=reshape(linspace(1,NPAR.nel,NPAR.nel),MSH.nz,nx);
ixs=[];
for iel=1:NPAR.nel
    [iz,ix]=find(e2ij==iel);
    fuel = Teff(iz,ix);
    mode = Rmod(iz,ix);
    [ixs]=InterpolateXS(XS.Data, MSH.z2mat(iel), fuel, mode, XS.Tf, XS.RhoMod, iel, ixs);
    if(nargin>5)
        tramp=XS.sav.tramp;
        if(time>=tramp)
            frod=0;
            funr=1;
        else
            frod=1-time/tramp;
            funr=time/tramp;
        end
        i1=find(loc(:,1)==iz);
        j1=find(loc(:,2)==ix);
        
        if(~isempty(i1) & ~isempty(j1))
%         if(iz==loc(1,1) & ix==loc(1,2))
%             disp('in loc');
%             pause
            ir = XS.imatR(iz,ix);
            iu = XS.imatU(iz,ix);
            ixs1=[];
            [ixs1]=InterpolateXS(XS.sav.Data, ir, fuel, mode, XS.Tf, XS.RhoMod, 1, ixs1);
            [ixs1]=InterpolateXS(XS.sav.Data, iu, fuel, mode, XS.Tf, XS.RhoMod, 2, ixs1);

            group = 1 ;
            ixs.rem(iel,group)   = ixs1.rem(1,group)*frod  + ixs1.rem(2,group)*funr ;
            ixs.diff(iel,group)  = ixs1.diff(1,group)*frod + ixs1.diff(2,group)*funr ;
            ixs.fiss(iel,group)  = ixs1.fiss(1,group)*frod + ixs1.fiss(2,group)*funr ;
            ixs.nuf(iel,group)   = ixs1.nuf(1,group)*frod  + ixs1.nuf(2,group)*funr ;
            ixs.invvelo(iel,group)  = ixs1.invvelo(1,group)*frod + ixs1.invvelo(2,group)*funr ;

            group = 2 ;
            ixs.rem(iel,group)   = ixs1.rem(1,group)*frod  + ixs1.rem(2,group)*funr ;
            ixs.diff(iel,group)  = ixs1.diff(1,group)*frod + ixs1.diff(2,group)*funr ;
            ixs.fiss(iel,group)  = ixs1.fiss(1,group)*frod + ixs1.fiss(2,group)*funr ;
            ixs.nuf(iel,group)   = ixs1.nuf(1,group)*frod  + ixs1.nuf(2,group)*funr ;
            ixs.invvelo(iel,group)  = ixs1.invvelo(1,group)*frod + ixs1.invvelo(2,group)*funr ;

            ixs.scatt(iel)  = ixs1.scatt(1)*frod + ixs1.scatt(2)*funr ;

        end
    end
end
