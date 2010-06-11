function XS=interpXS(time,XS);


tramp=XS.sav.tramp;

if(time>=tramp)
    frod=0;
    funr=1;
else
    frod=1-time/tramp;
    funr=time/tramp;
end

XS.dif(4,:) = frod * XS.sav.dif(4,:) + funr * XS.sav.dif(1,:) ;

XS.rem(4,:) = frod * XS.sav.rem(4,:) + funr * XS.sav.rem(1,:) ;

XS.nuf(4,:) = frod * XS.sav.nuf(4,:) + funr * XS.sav.nuf(1,:) ;

XS.sca(4,:) = frod * XS.sav.sca(4,:) + funr * XS.sav.sca(1,:) ;

% chi, b, btot, invv, lambda, nu, nup, nud : unchanged

XS.fis(4,:) =XS.nuf(4,:)./XS.nu(4,:);

XS.nudf(4,:)=XS.nud(4,:) .* XS.fis(4,:) ;
XS.nupf(4,:)=XS.nup(4,:) .* XS.fis(4,:) ; 