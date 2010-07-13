function [ixs]=InterpolateXS(Data, matnum, fuel, mode, T_f, RhoMod, iel, ixs)

% T_f = [.5000000E+03  .7602200E+03  .8672700E+03  .9218800E+03  .1500000E+04] ;
% RhoMod = [.6413994E+03  .7114275E+03  .7694675E+03  .7724436E+03  .7813064E+03 .8100986E+03] ;

fuel = fuel + 273.15 ;

% safeguards
eps=0.1;
fuel=min(fuel, T_f(end)-eps);
fuel=max(fuel, T_f(1)+eps  );
mode=min(mode, RhoMod(end)-eps);
mode=max(mode, RhoMod(1)+eps  );

[j,k] = GetIndices(fuel, mode, T_f, RhoMod) ;

group = 1 ;
ixs.rem(iel,group)  = InterpolateData(Data{matnum}.rem{group} , T_f, RhoMod, fuel, mode, j, k) ;
ixs.diff(iel,group) = InterpolateData(Data{matnum}.diff{group}, T_f, RhoMod, fuel, mode, j, k) ;
ixs.fiss(iel,group) = InterpolateData(Data{matnum}.fiss{group}, T_f, RhoMod, fuel, mode, j, k) ;
ixs.nuf(iel,group)  = InterpolateData(Data{matnum}.nuf{group} , T_f, RhoMod, fuel, mode, j, k) ;
ixs.invvelo(iel,group) = Data{matnum}.InvVelo(group);

group = 2 ;
ixs.rem(iel,group)  = InterpolateData(Data{matnum}.rem{group} , T_f, RhoMod, fuel, mode, j, k) ;
ixs.diff(iel,group) = InterpolateData(Data{matnum}.diff{group}, T_f, RhoMod, fuel, mode, j, k) ;
ixs.fiss(iel,group) = InterpolateData(Data{matnum}.fiss{group}, T_f, RhoMod, fuel, mode, j, k) ;
ixs.nuf(iel,group)  = InterpolateData(Data{matnum}.nuf{group} , T_f, RhoMod, fuel, mode, j, k) ;
ixs.invvelo(iel,group) = Data{matnum}.InvVelo(group);

ixs.scatt(iel) = InterpolateData(Data{matnum}.scatt, T_f, RhoMod, fuel, mode, j, k) ;


return ;