function myXS =load_XS_all( read_entire_file, z2med)


if(read_entire_file)
    fid = fopen('nemtab');
    %     C = textscan(fid,'%f', 'delimiter', '  ', ...
    %         'headerLines', 6, 'commentStyle', '*') ;
    C = textscan(fid,'%f', 'delimiter', '  ', 'emptyValue', 0, ...
        'headerLines', 6, 'commentStyle', '*') ;
    fclose(fid) ;
    Data = C{1} ;
    clear C
    index = find(Data) ;
    tfnum = 5 ;  rhomodnum = 6 ;
    nmat_unrodded = 438;
    sign = +1;
    [XS.Data,XS.Tf,XS.RhoMod] = ReadXSData(index, Data, tfnum, rhomodnum, nmat_unrodded, sign) ;

    [a,b] = find(z2med==0);
    if(~isempty(a)) % not empty
        error('no materials should have a zero #')
    end

    my_ind  = 0;
    my_imat = z2med*0;

    [n1 n2]=size(z2med);
    my_imatU = zeros(n1,n2);
    my_imatR = zeros(n1,n2);
    
    z2medORI=z2med;
    z2med=abs(z2med);
    z=reshape(z2med,n1*n2,1);
    
    im2 = max(z);   
    % there are some unrodded compositions
    if(im2>0) 
        im1 = min(z);
        ind = linspace(im1,im2,im2-im1+1);
        for i=1:length(ind)
            [a,b] = find(z==ind(i));
            if(~isempty(a)) % not empty
                my_ind = my_ind + 1;
                myXS.Data{my_ind} = XS.Data{ind(i)};
                [aa,bb] = find(z2medORI==ind(i));
                if( length(a)==length(aa) & length(b)==length(bb) )
                    if(a==aa & b==bb)
                        for j1=1:length(a)
                            my_imat(a(j1),b(j1)) = my_ind;
                            my_imatU(a(j1),b(j1)) = my_ind;
                        end
                    else
                        error('a ~= aa & b~=bb')
                    end
                end
                % the unrodded version of a rodded material
                [aa,bb] = find(z2medORI==-ind(i));
                if( length(a)==length(aa) & length(b)==length(bb) )
                    if(a==aa & b==bb)
                        for j1=1:length(a)
                            my_imatU(a(j1),b(j1)) = my_ind;
                        end
                    else
                        error('a ~= aa & b~=bb')
                    end
                end
            end
        end
        myXS.Tf     = XS.Tf;
        myXS.RhoMod = XS.RhoMod;
    end

    z2med=-z2med; 
    z=-z;

    % rodded compositions
    im1 = min(z);
    if(im1<0)
        fid = fopen('nemtabr_3.01');
        C = textscan(fid,'%f', 'delimiter', '  ', 'emptyValue', 0, ...
            'headerLines', 6, 'commentStyle', '*') ;
        fclose(fid) ;
        DataR = C{1} ;
        clear C ;
        index = find(DataR) ;
        tfnum = 5 ;  rhomodnum = 6 ;
        nmat_rodded = 195;
        sign = -1;
        [XSR.Data,XS.Tf,XS.RhoMod] = ReadXSData(index, DataR, tfnum, rhomodnum, nmat_rodded, sign) ;

        im2 = max(z);
        ind = linspace(im1,im2,im2-im1+1);
        for i=1:length(ind)
            [a,b] = find(z2med==ind(i));
            if(~isempty(a)) % not empty
                my_ind = my_ind + 1;
                myXS.Data{my_ind} = XSR.Data{-ind(i)}; 
                [aa,bb] = find(z2medORI==ind(i));
                if( length(a)==length(aa) & length(b)==length(bb) )
                    if(a==aa & b==bb)
                        for j1=1:length(a)
                            my_imat(a(j1),b(j1)) = my_ind;
                            my_imatR(a(j1),b(j1)) = my_ind;
                        end
                    else
                        error('a ~= aa & b~=bb')
                    end
                end
                % unrodded version of a currently rodded composition
                [aa,bb] = find(z2medORI==-ind(i));
                if( length(a)==length(aa) & length(b)==length(bb) )
                    if(a==aa & b==bb)
                        for j1=1:length(a)
                            my_imatR(a(j1),b(j1)) = my_ind;
                        end
                    else
                        error('a ~= aa & b~=bb')
                    end
                end
%                 my_ind = my_ind + 1;
%                 myXS.Data{my_ind} = XS.Data{-ind(i)};
            end
        end
        myXS.Tf     = XS.Tf;
        myXS.RhoMod = XS.RhoMod;
        clear XSR ;
    end
    clear XS ;
    
    myXS.nmat = length(myXS.Data);
    myXS.imat = my_imat;
    myXS.imatU = my_imatU;
    myXS.imatR = my_imatR;
    myXS.z2medORI=z2medORI;
    save my_XS.mat myXS ;
    %     [my_imat my_imatU my_imatR z2medORI]
else
    load my_XS.mat ;
end


% create the chi array
myXS.chi=zeros( myXS.nmat ,2 );
myXS.chi(:,1) = 1;

% # of energy groups
myXS.G = length(myXS.chi(1,:));

% kin lambda
myXS.lambda        = zeros(4,6);
myXS.lambda(1,1:6) = [0.0124 0.0305 0.111 0.301 1.14 30.1];

% kin beta
myXS.b        = zeros(4,6);
myXS.b(1,1:6) = [21 142 127 257 75 27]*1e-5;

for i=2:myXS.nmat
    myXS.lambda(i,1:6)= myXS.lambda(1,1:6);
    myXS.b(i,1:6)     = myXS.b(1,1:6);
end

% beta tot
for i=1:myXS.nmat
    myXS.btot(i)     = sum(myXS.b(i,1:6));
end

% # of precursor groups
bidon = size(myXS.lambda);
myXS.nprec=bidon(2);

%
myXS.tramp = .250;
% myXS.tramp = .010;

myXS.sav  = myXS;

