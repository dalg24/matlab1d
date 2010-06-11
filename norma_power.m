function [Pow,P0]=norma_power(Pow,Watt0,impr)

factor = sum(Pow);

if(abs(factor)>1e-12)
    Pow = Pow /factor * Watt0;
end

P0 = sum(Pow) ;

[Pmax,i]=max(Pow);

if(impr)
    disp(sprintf('Max power %d (W) at axial node %i ',Pmax,i));
end
