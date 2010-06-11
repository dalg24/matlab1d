function [out]=InterpolateData(c, Tf, RhoMod, fuel, mod, j, k)

tflen = length(Tf) ;
rholen = length(RhoMod) ;

if (j < tflen && k < rholen)
    
%     if(fuel<Tf(j)),  fuel=Tf(j)  ,end
%     if(fuel>Tf(j+1)),fuel=Tf(j+1),end
%     if(mod<RhoMod(k)),  mod=RhoMod(k)  ,end
%     if(mod>RhoMod(k+1)),mod=RhoMod(k+1),end
    % fuel and mod within limits
    t = (fuel - Tf(j))/(Tf(j+1) - Tf(j)) ;
    u = (mod - RhoMod(k))/(RhoMod(k+1) - RhoMod(k)) ;

    out = (1 - t)*(1 - u)*c(k,j) + t*(1 - u)*c(k,j+1) + t*u*c(k+1,j+1) + (1 - t)*u*c(k+1,j) ;
    
else
    error('Out of limits') ;
    
    if(j >= tflen && k >= rholen)

        out = c(rholen,tflen) ;

    elseif(j >= tflen)

        u = (mod - RhoMod(k))/(RhoMod(k+1) - RhoMod(k)) ;

        out = (1 - u)*c(k,tflen) + u*c(k+1,tflen) ;

    else

        t = (fuel - Tf(j))/(Tf(j+1) - Tf(j)) ;

        out = (1 - t)*c(rholen,j) + t*c(rholen,j+1) ;

    end
end