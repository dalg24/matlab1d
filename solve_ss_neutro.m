function [keff,v]=solve_ss_neutro(MSH,NPAR,XS,BC,therm,hydro,vec)

[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);
[P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod);
%     [kinf]=comp_kinf(ixs); figure(12); plot(kinf);hold on


G=XS.G;
% call matlab's eigensolver
OPTS.disp=0;
[v,d]=eigs(inv(M)*P,1,'lm',OPTS);
% make sure the fundamental mode is positive
a=length(find(v>0));
b=length(find(v<0));
if(b>a)
    v(:,1)=-v(:,1);
end

% keff
keff=max(diag(d));


% % % whether the fission XS should be normalized by keff
% % if(log_renorm)
% %     %     for i=1:length(XS.Data)
% %     %         for g=1:G
% %     %             XS.Data{i}.nuf{g}(:,:)=XS.Data{i}.nuf{g}(:,:)/keff;
% %     %         end
% %     %     end
% % 
% % 
% %     for ii=1:length(XS.Data)
% %         if(length(XS.Data{ii})~=0)
% %             for g=1:XS.G
% %                 XS.Data{ii}.nuf{g}(:,:)=XS.Data{ii}.nuf{g}(:,:)/keff;
% % %                 XS.sav.Data{ii}.nuf{g}(:,:)=XS.sav.Data{ii}.nuf{g}(:,:)/keff;
% %             end
% %         end
% %     end
% % 
% %     XS.sav=[];
% %     XS.sav=XS;
% % end


