function [xn,IV,A,NPAR]=init_kin(XS,NPAR,MSH,BC,therm,hydro,vec)

% init kin
n = NPAR.dof;
G = XS.G;
L = XS.nprec;

prec_dof = L*n;
NPAR.prec = NPAR.neu  + prec_dof;
NPAR.ther = NPAR.ther + prec_dof;
NPAR.siz  = NPAR.siz  + prec_dof;
xn = zeros(NPAR.siz,1);

% copy init data into beginning-of-time-step solution vector
xn(1:NPAR.neu)            = vec(1:NPAR.neu);
xn(NPAR.prec+1:NPAR.ther) = vec(NPAR.neu+1:NPAR.ther-prec_dof);
xn(NPAR.ther+1:NPAR.siz)  = vec(NPAR.ther-prec_dof+1:NPAR.siz-prec_dof);

[Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xn );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xn);

% create 1/V
[IV,ixs]=create_mat_invv(XS,MSH,NPAR,Teff,Rmod);

time=0.;
[A,ixs]=create_mat_tn(MSH,XS,NPAR,BC,Teff,Rmod,time,IV);


B =A( G*n+1:(G+L)*n, 1:G*n );
La=A( G*n+1:(G+L)*n, G*n+1:(G+L)*n);

Cn = -inv(La)*B*xn(1:NPAR.neu);
xn(NPAR.neu+1:NPAR.prec) = Cn;

iplot=0;
if(iplot>0)
    figure(2);
    leng=length(z2med)*dx;
    xx=(linspace(0,leng,n))';
    for i=1:L,
        plot(xx,Cn((i-1)*n+1:i*n))
        hold on
    end
end

if(iplot>10)
    % verif
    X=A( 1:G*n , G*n+1:(G+L)*n  );
    Pp_M=A(1:G*n,1:G*n);
    qq=Pp_M*v+X*Cn;
    figure(3)
    plot(qq)

    figure(4);
    plot(xx,v(1:n,1)/max(v(1:n,1)) ,xx,v(n+1:2*n,1)/max(v(n+1:2*n,1)));
    % figure(5);
    % plot(xx,v(1:n,1)/max(v(1:n,1)) ,xx,v(n+1:2*n,1)/max(v(n+1:2*n,1)));
    % ,...
    figure(5);
    plot(xx,Cn(0*n+1:1*n,1)/max( Cn(0*n+1:1*n,1) ),...
        xx,Cn(1*n+1:2*n,1)/max( Cn(1*n+1:2*n,1) ),...
        xx,Cn(2*n+1:3*n,1)/max( Cn(2*n+1:3*n,1) ),...
        xx,Cn(3*n+1:4*n,1)/max( Cn(3*n+1:4*n,1) ),...
        xx,Cn(4*n+1:5*n,1)/max( Cn(4*n+1:5*n,1) ),...
        xx,Cn(5*n+1:6*n,1)/max( Cn(5*n+1:6*n,1) ) )

    figure(6);
    plot(xx,v(1:n,1)/max(v(1:n,1)) ,xx,v(n+1:2*n,1)/max(v(n+1:2*n,1)),...
        xx,Cn(0*n+1:1*n,1)/max( Cn(0*n+1:1*n,1) ),...
        xx,Cn(1*n+1:2*n,1)/max( Cn(1*n+1:2*n,1) ),...
        xx,Cn(2*n+1:3*n,1)/max( Cn(2*n+1:3*n,1) ),...
        xx,Cn(3*n+1:4*n,1)/max( Cn(3*n+1:4*n,1) ),...
        xx,Cn(4*n+1:5*n,1)/max( Cn(4*n+1:5*n,1) ),...
        xx,Cn(5*n+1:6*n,1)/max( Cn(5*n+1:6*n,1) ) )
    legend('f','t','1','2','3','4','5','6')
    % Cbis=-inv(X)*Pp_M*v;
    % for i=1:L,
    %     figure(3+i);
    %     plot(xx,Cn((i-1)*n+1:i*n) - Cbis((i-1)*n+1:i*n) )
    %     hold on
    % end
end