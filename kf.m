function [out]=kf(t,t0);
% C=======================================================================
% C
% C     fuel conductivity evaluation (or conductivity integral evalutation)
% C     units: W/(degre.m)
% C
% C=======================================================================
% *>t   Temperature in degres Celsius
% *>t0  Temperature in degres Celsius
% C=======================================================================
% C
% C if t0 is present, the conductivity integral is computed between the
% C bounds t and to: i.e. integral (from t0 to t) of k(s)ds
% C if only t is provided, k(t) is computed 
% C
% C=======================================================================

c(1,1)=   40.4000e+0;
c(2,1)=  464.0000e+0;
c(3,1)=    1.2160e-4;
c(4,1)=    1.8670e-3;
c(5,1)=    0.0191e+0;

c(1,2)=   33.0000e+0;
c(2,2)=  375.0000e+0;
c(3,2)=    1.5400e-4;
c(4,2)=    1.7100e-3;
c(5,2)=    0.0171e+0;

TLIM(1)= 1650.e+0;
TLIM(2)= 1550.e+0;
% C
% C=======================================================================
% C
COMPO=1;
FD=0.95;

if(COMPO==1)
    BETA=2.58e+0 -0.58e-3*t;
    P=(1.-BETA*(1-FD))/(1.-BETA*0.05);
else
    P=1.10125e+0*FD/(2.43e+0-1.43e+0*FD);
end

if(nargin==1)
    if (t<TLIM(COMPO)) 
        k=100*P*( c(1,COMPO)/(c(2,COMPO)+t) + c(3,COMPO)*exp(c(4,COMPO)*t) );
    else
        k=100*P*( c(5,COMPO)+ c(3,COMPO)*exp(c(4,COMPO)*t) );
    end
elseif(nargin==2)
    if (t<TLIM(COMPO)) 
        int=100*P*( c(1,COMPO)*log((c(2,COMPO)+t)/(c(2,COMPO)+t0)) +...
            c(3,COMPO)/c(4,COMPO)*(exp(c(4,COMPO)*t)-exp(c(4,COMPO)*t0)) );
    else
        int=100*P*( c(5,COMPO)*(t-t0)+ ...
            c(3,COMPO)/c(4,COMPO)*( exp(c(4,COMPO)*t) - exp(c(4,COMPO)*t0) ) );
    end
end

if(nargin==1)
    out=k;
else
    out=int;
end

% C=======================================================================
% *>TEMP   Temperature locale en degres Celsius
% C>COMPO  UO2 ou MOX (1 ou 2)
% C>FD     fraction de densite theorique (typiquement 0.95)
% C=======================================================================

