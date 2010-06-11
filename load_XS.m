function myXS =load_XS( read_entire_file, z2med)


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
    z=reshape(z2med,n1*n2,1);
    
%     im1 = min(min(z2med));
%     im2 = max(max(z2med));
    im1 = min(z);
    im2 = max(z);
    if(im2>0)
        % im1 = min(min(abs(z2med)));
        qq=find(z>0);
        im1 = min(z(qq));
        ind = linspace(im1,im2,im2-im1+1);
        for i=1:length(ind)
            [a,b] = find(z2med==ind(i));
            if(~isempty(a)) % not empty
                my_ind = my_ind + 1;
                myXS.Data{my_ind} = XS.Data{ind(i)};
                if(length(a)~=length(b))
                    error('len a ~= len b')
                end
                for j1=1:length(a)
                    my_imat(a(j1),b(j1)) = my_ind;
                end
            end
        end
        myXS.Tf     = XS.Tf;
        myXS.RhoMod = XS.RhoMod;
%         clear XS ;
    end

%     im1 = min(min(z2med));
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

        qq=find(z<0);
        im2 = max(z(qq));
%         im2 = -max(max(-z2med));
        ind = linspace(im1,im2,im2-im1+1);
        for i=1:length(ind)
            [a,b] = find(z2med==ind(i));
            if(~isempty(a)) % not empty
                my_ind = my_ind + 1;
                myXS.Data{my_ind} = XSR.Data{-ind(i)};
                if(length(a)~=length(b))
                    error('len a ~= len b')
                end
                for j1=1:length(a)
                    my_imat(a(j1),b(j1)) =my_ind;
                end
                my_ind = my_ind + 1;
                myXS.Data{my_ind} = XS.Data{-ind(i)};
            end
        end
        myXS.Tf     = XS.Tf;
        myXS.RhoMod = XS.RhoMod;
%         clear XS ;
    end
    clear XS ;
    
    myXS.nmat = length(myXS.Data);
    myXS.imat = my_imat;
    save my_XS.mat myXS ;
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
myXS.tramp = .00000000000000000000250;
myXS.tramp = .250;

myXS.sav  = myXS;

