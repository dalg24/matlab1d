function [XS,Tf,RhoMod]=ReadXSData(index, Data, tfnum, rhomodnum, nmat, sign)

Tf = [.5000000E+03  .7602200E+03  .8672700E+03  .9218800E+03  .1500000E+04] ;
RhoMod = [.6413994E+03  .7114275E+03  .7694675E+03  .7724436E+03  .7813064E+03 .8100986E+03] ;

skip    = tfnum + rhomodnum ;
readnum = tfnum * rhomodnum ;

% be = skip ;    en = skip+readnum ;

for mm=1:nmat,
    if(mm==1)
        mat = Data(index(1)) ;
        i = skip + 2 ;
    else
        i = i+3 ;
        mat = Data(index(i-1)) ;
        i = i + skip ;
    end
    %     if(mm~=mat)
    %         error('mm~=mat')
    %     end

    % group #1
    %%%%%%%%%%
    group = 1 ;
    for row = 1 : rhomodnum
        for col = 1 : tfnum
            XS{mat}.diff{group}(row,col) = Data(index(i)) ;
            i = i + 1 ;
        end
    end

    i = i + skip ;
    for row = 1 : rhomodnum
        for col = 1 : tfnum
            XS{mat}.rem{group}(row,col) = Data(index(i)) ;
            i = i + 1 ;
        end
    end

    if((sign>0)&(mat==24|mat==25|mat==26))
        i = i + skip ;
        XS{mat}.fiss{group}(1 : rhomodnum,1 : tfnum) = 0;
        i = i + skip ;
        XS{mat}.nuf{group}(1 : rhomodnum,1 : tfnum)  = 0;
    else
        i = i + skip ;
        for row = 1 : rhomodnum
            for col = 1 : tfnum
                XS{mat}.fiss{group}(row,col) = Data(index(i)) ;
                i = i + 1 ;
            end
        end

        i = i + skip ;
        for row = 1 : rhomodnum
            for col = 1 : tfnum
                XS{mat}.nuf{group}(row,col) = Data(index(i)) ;
                i = i + 1 ;
            end
        end
    end

    i = i + skip ;
    for row = 1 : rhomodnum
        for col = 1 : tfnum
            XS{mat}.scatt(row,col) = Data(index(i)) ;
            i = i + 1 ;
        end
    end

    % complete the calculattion of rem xs
    XS{mat}.rem{group}(1:rhomodnum,1:tfnum) = XS{mat}.rem{group}(1:rhomodnum,1:tfnum) + ...
        XS{mat}.scatt(1:rhomodnum,1:tfnum);

    % group #2
    %%%%%%%%%%
    group = 2 ;

    i = i + skip ;
    for row = 1 : rhomodnum
        for col = 1 : tfnum
            XS{mat}.diff{group}(row,col) = Data(index(i)) ;
            i = i + 1 ;
        end
    end

    i = i + skip ;
    for row = 1 : rhomodnum
        for col = 1 : tfnum
            XS{mat}.rem{group}(row,col) = Data(index(i)) ;
            i = i + 1 ;
        end
    end

    if((sign>0)&(mat==24|mat==25|mat==26))
        i = i + skip ;
        XS{mat}.fiss{group}(1 : rhomodnum,1 : tfnum) = 0;
        i = i + skip ;
        XS{mat}.nuf{group}(1 : rhomodnum,1 : tfnum)  = 0;
    else
        i = i + skip ;
        for row = 1 : rhomodnum
            for col = 1 : tfnum
                XS{mat}.fiss{group}(row,col) = Data(index(i)) ;
                i = i + 1 ;
            end
        end

        i = i + skip ;
        for row = 1 : rhomodnum
            for col = 1 : tfnum
                XS{mat}.nuf{group}(row,col) = Data(index(i)) ;
                i = i + 1 ;
            end
        end
    end

    % skipping Xe always
    if((sign>0)&(mat==24|mat==25|mat==26))
        i = i + skip + 0 ;
    else
        i = i + skip + readnum ;
    end
    XS{mat}.InvVelo = [Data(index(i)) Data(index(i+1))] ;
end
% %     % Rodded material
% %     i = i+3 ;
% %     mat = Data(index(i-1)) ;
% %
% %     % group #1
% %     %%%%%%%%%%
% %     group = 1 ;
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.diff{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.rem{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.fiss{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.nuf{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.scatt(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     % group #2
% %     %%%%%%%%%%
% %     group = 2 ;
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.diff{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.rem{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.fiss{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip ;
% %     for row = 1 : rhomodnum
% %         for col = 1 : tfnum
% %             XS{mat}.nuf{group}(row,col) = Data(index(i)) ;
% %             i = i + 1 ;
% %         end
% %     end
% %
% %     i = i + skip + readnum ;
% %     XS{mat}.InvVelo = [Data(index(i)) Data(index(i+1))] ;