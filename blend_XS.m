function XS = blend_XS( XS, MSH, w )

% [n1,n2] = size(loc);
% if(n2~=MSH.ndim)
%     error('wrong size loc')
% end

[a,b]=find(MSH.z2med<0);

for k=1:length(a),
    
    % material number
    ir = XS.imatR(a(k),b(k));
    iu = XS.imatU(a(k),b(k));
    
    % XS.Data{ir} = XS.sav.Data{ir} * w +  XS.sav.Data{iu} * (1-w);     
    XS.Data{ir}.diff{1,1}(:,:) = XS.sav.Data{ir}.diff{1,1}(:,:) * w +  XS.sav.Data{iu}.diff{1,1}(:,:) * (1-w);
    XS.Data{ir}.diff{1,2}(:,:) = XS.sav.Data{ir}.diff{1,2}(:,:) * w +  XS.sav.Data{iu}.diff{1,2}(:,:) * (1-w);
    XS.Data{ir}.rem{1,1}(:,:)  = XS.sav.Data{ir}.rem{1,1}(:,:)  * w +  XS.sav.Data{iu}.rem{1,1}(:,:)  * (1-w);
    XS.Data{ir}.rem{1,2}(:,:)  = XS.sav.Data{ir}.rem{1,2}(:,:)  * w +  XS.sav.Data{iu}.rem{1,2}(:,:)  * (1-w);
    XS.Data{ir}.fiss{1,1}(:,:) = XS.sav.Data{ir}.fiss{1,1}(:,:) * w +  XS.sav.Data{iu}.fiss{1,1}(:,:) * (1-w);
    XS.Data{ir}.fiss{1,2}(:,:) = XS.sav.Data{ir}.fiss{1,2}(:,:) * w +  XS.sav.Data{iu}.fiss{1,2}(:,:) * (1-w);
    XS.Data{ir}.nuf{1,1}(:,:)  = XS.sav.Data{ir}.nuf{1,1}(:,:)  * w +  XS.sav.Data{iu}.nuf{1,1}(:,:)  * (1-w);
    XS.Data{ir}.nuf{1,2}(:,:)  = XS.sav.Data{ir}.nuf{1,2}(:,:)  * w +  XS.sav.Data{iu}.nuf{1,2}(:,:)  * (1-w);

    XS.Data{ir}.scatt(:,:)     = XS.sav.Data{ir}.scatt(:,:) * w +  XS.sav.Data{iu}.scatt(:,:) * (1-w);


    XS.Data{ir}.InvVelo(:,:) = XS.sav.Data{ir}.InvVelo(:,:) * w +  XS.sav.Data{iu}.InvVelo(:,:) * (1-w);
    XS.Data{ir}.InvVelo(:,:) = XS.sav.Data{ir}.InvVelo(:,:) * w +  XS.sav.Data{iu}.InvVelo(:,:) * (1-w);
    
end

    

XS.sav = XS ;
