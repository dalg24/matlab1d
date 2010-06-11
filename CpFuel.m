function [cp]= CpFuel(T)
	TK=T+273.15E+0;
	cp = 162.3E+0 +0.3038E+0*TK -2.391E-4*TK^2 +6.404E-8*TK^3;
end % functio CpFuel

