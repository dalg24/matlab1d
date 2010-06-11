function [j, k]=GetIndices(fuel, mode, T_f, RhoMod)

j = 1 ; k = 1 ;
index = find(T_f<=fuel) ;
if index == 0
    j = 1 ;
else
    j = index(end) ;
end
index = find(RhoMod<=mode) ;
if index == 0
    k = 1 ;
else
    k = index(end) ;
end