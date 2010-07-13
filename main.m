clc; close all; %clear all;    

%%%%%%%%% input %%%%%%%%%%
% axial refinement
ref_z=1;
% load input data
[MSH, XS, therm, hydro] = loadinput( ref_z );

%%%%%%%%% numerical parameters %%%%%%%%%%%
[NPAR, BC] = numpar( MSH, XS.G );

%%%%%%%%% initialize solution vector %%%%%
[vec] = init_vec( MSH, XS.G, hydro.T_in, hydro.h_in, NPAR );

% initial guess for jeff
keff=1;

% if(MSH.dim==1)
%     load save1r.mat;
% end

% figure ID for plotting
fid=99;
tol=1e-6;

if(NPAR.NLNewton)

    % compute the initial steady-state fission source 
    [fisrc_old,P,M,ixs]=compute_ss_fiss_src( MSH, NPAR, BC, XS, therm, hydro, vec);

    for k=1:NPAR.max_ss_iter,
        
        % plot initial Teff and Tmod
        my_TH_plot(fid, true, true, MSH, NPAR, therm, hydro, vec );

        disp(sprintf('Steady-state N-TH iteration # %4i',k));

        % compute residula (the fission src term is fixed)
        [res,err2] = comp_residual( MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old );
        if(NPAR.MatFree)
            % restart=10; % gmres restart
            tol = min(1e-1,1e-2*err2); tol = max(tol, 1e-6);   % gmres tol
            maxit = 80;
            epsi  = 1e-7;
        else
            % compute the analytical Jacobian
            J=compute_num_jacobian( res, MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old);
        end

        % perform the Newton solve for the N-TH system with a fixed fission source
        for i=1:NPAR.MaxNewton
            if(NPAR.MatFree)
                [delta_vec,flag,relres,iter] = gmres( @comp_res_diff, -res, [],tol,maxit,[],[],[], ...
                    MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old, res, epsi);
                if(NPAR.PreCond)
                    [delta_vec,bidon]=OS_ss(MSH, NPAR, XS, therm, hydro, BC, delta_vec, vec);
                end
            else
                delta_vec =-J\res;
            end
            vec = vec + delta_vec;
            [res,err2] = comp_residual( MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old );
            % compute errors
            % err=norm(delta_vec,inf);
            ii=find(vec>0);
            dv_v = norm( delta_vec(ii)./vec(ii), inf) ;
            % printouts
            disp(sprintf('      Newton iteration %i,  max(delta/vector) %g, max(res/vector) %g, tol %g',i,dv_v,err2,tol));
            %             if( err2 <1e-4 | i==NPAR.MaxNewton ), pause, end
            if( err2 <tol)
                % plot Teff and Tmod
                my_TH_plot(fid+1, true, false, MSH, NPAR, therm, hydro, vec );
                disp('      *** converged nonlinear solve ***');disp(' ')
                break
            else
                if(NPAR.MatFree)
                    [delta_vec,flag,relres,iter] = gmres( @comp_res_diff,-res,[],tol,maxit,[],[],[], MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old, res, epsi);
                    if(NPAR.PreCond)
                        [delta_vec,bidon]=OS_prec_ss(MSH, NPAR, XS, therm, hydro, BC, delta_vec, vec);
                    end
                else
                    % compute the analytical Jacobian
                    J=compute_num_jacobian( res, MSH, XS, NPAR, BC, therm, hydro, vec, fisrc_old);
                end
                % plot Teff and Tmod
                my_TH_plot(fid, true, true, MSH, NPAR, therm, hydro, vec );
            end
        end

        % solve ss neutronics to get keff with the TH conditions determined
        % from the above Newton solve
        keff_old = keff;
        [keff,vec(1:NPAR.neu)] = solve_ss_neutro(MSH,NPAR,XS,BC,therm,hydro,vec);

        % update the steady-state fission source
        [fisrc_old,P,M,ixs]=compute_ss_fiss_src( MSH, NPAR, BC, XS, therm, hydro, vec);

        % compute power
        [Pow,Plevel(1)] = comp_power( ixs, MSH, NPAR, vec(1:NPAR.neu), XS.G );
        [Pow,Plevel(1)] = norma_power( Pow, therm.Watt0, false );
%         Plevel(1)

        % error in the eigenvalue
        err_eig = abs(keff-keff_old)*1e5;
        disp(sprintf('      Steady-state coupled outer , iter: %i, eigenvalue estimate: %6d, err_eig %d',...
            k,keff,err_eig));
        disp(' ');disp(' ');
        if(err_eig<1e-2)
            disp('* * * * converged eigenproblem using Newton solve * * * *');
            break
        end

    end
    
else
    
    for i=1:NPAR.max_ss_iter

        disp(sprintf('OS iteration %i',i));
        
        [vec_new,keff,Plevel]=OS_ss(MSH, NPAR, XS, therm, hydro, BC, vec);
        
        error=norm(vec_new-vec,2);
        disp(sprintf('    error OS = %g',error));

        vec=NPAR.damp*vec_new+(1-NPAR.damp)*vec;       
        
        if(error<1e-6)
            disp('converged OS');
            break
        end
    end
end


% make critical
for ii=1:length(XS.Data)
    if(length(XS.Data{ii})~=0)
        for g=1:XS.G
            XS.Data{ii}.nuf{g}(:,:)=XS.Data{ii}.nuf{g}(:,:)/keff;
            XS.sav.Data{ii}.nuf{g}(:,:)=XS.sav.Data{ii}.nuf{g}(:,:)/keff;
        end
    end
end
% verify criticality
[Teff]       = extract_Teff( MSH, NPAR.neu, therm.rw, vec );
[Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, vec);
[P,M,ixs]    = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod);
fl=vec(1:NPAR.neu);
a=M*fl-P*fl;
[q,qq]=eigs(M\P);
max(max(qq))
clear q qq a;

my_TH_plot(fid+1, true, false, MSH, NPAR, therm, hydro, vec );


% [fl]=re_ordering(MSH,vec(1:NPAR.neu),NPAR.gn,XS.G,NPAR.dof);

if(MSH.dim==1)
%     save save1rjfnk.mat vec keff;
    save save_os_ss.mat vec keff;
% else
%     save save2r.mat vec keff;
end

error('eeeee')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NPAR.MatFree=false;
NPAR.PreCond=true;

tend=5; dt=0.01; nt=tend/dt; t=0.;
aux=abs(nt*dt-tend);
if(aux>1e-14)
    error('wrong nt')
end

% retrieve Plevel(1)
[Pow,Plevel(1)] = comp_power( ixs, MSH, NPAR, vec(1:NPAR.neu), XS.G );
% normalize fluxes
vec(1:NPAR.neu)= vec(1:NPAR.neu)/Plevel(1);
iplot=0;
% compute SS vector
XS.loc=[18,1]*0;
[xn,IV,An,NPAR]=init_kin(XS,NPAR,MSH,BC,therm,hydro,vec);
xnp1 = xn;
clear vec;

% compute power with normalized vector
[Pow,Plevel(1)] = comp_power( ixs, MSH, NPAR, xn(1:NPAR.neu), XS.G );
[Pow,Plevel(1)] = norma_power( Pow, therm.Watt0, true );

% save old values
Pown=Pow;
[Rmodn, Tmodn] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xn);
[flxthn] = extract_flxth( MSH, hydro, xn(NPAR.prec+1:NPAR.ther), Tmodn, ...
    xn(NPAR.ther+1:NPAR.siz) );
% % flxthn(end-1:end)=0;

ns=length(diag(IV));

% filename = input('filename for transient = ');
% test 1 : forgot to have xnp1 to compute anp1, hence old temop all the time
% test 2 : corrected previous mistake
% test 3 : only fuel+neutro
% test 4 : added update of flxth after hydro
% test 4u: OS with only 1 iter
% test 5u: OS with only 1 iter, neutro last operator action


% need to run a new OS test b/c of hc/hcn in solve_tf_hc
% new geom conf
name = 'Case1_MF_prec_tol4_5sec_cv_';
name = 'Case1_OS_100it_5sec_cv_';
name = 'test';
time_str=sprintf('%0.4f',dt);
time_str(1:2)=[];
ncut=4-ceil(-log10(dt));
time_str(4-ncut+1:4)=[];
filename=[name(1:end),time_str,'.txt'];
disp(filename)


fid=fopen(filename,'w');
aux=[dt,tend,nt]; fwrite(fid,aux,'double'); clear aux;
fwrite(fid,xnp1,'double');

save structures.mat MSH NPAR therm hydro BC XS;

XS.loc=[20,1;21,1;22,1]*1;

% XS.loc=[39,1; 40,1; 41,1; 42,1; 43,1; 44,1; 45,1; 46,1; 47,1; 48,1]*1;
do.neutro=1;
do.thermo=1;
do.hydro =1;

test.neutro=0;
test.thermo=0;
test.hydro =0;
flxth0=flxthn;

t_ramp=min(1+linspace(0,tend,nt+1),2);
NPAR.method='CN';


neutro_first=true;

therm.Plevel = Plevel;
if(NPAR.MatFree)
    [resn,err2] = comp_residual_tr( MSH, XS, NPAR, BC, therm, hydro, xn, t, dt, IV );
    max_os_iter = 0;
else
    max_os_iter = 100;
end

iter_=zeros(1,2);
for it=1:nt,
    disp(sprintf('n-D transient: time=%d, tend=%d',t,tend))



    if(NPAR.MatFree)

        [resn  ,err2] = comp_residual_tr( MSH, XS, NPAR, BC, therm, hydro, xn  , t,    dt, IV      );
        [resnp1,err2] = comp_residual_tr( MSH, XS, NPAR, BC, therm, hydro, xnp1, t+dt, dt, IV, resn);

        for i=1:NPAR.MaxNewton
            %         restart=10; % gmres restart
%             tol=min(1e-1,1e-2*err2);
%             tol=max(tol, 1e-6);   % gmres tol
            tol=1e-4;
            maxit=NPAR.siz;
            epsi = 1e-7;
            [delta_vec,flag,relres,iter] = gmres( @comp_res_diff_tr,-resnp1,[],tol,maxit,[],[],[], ...
                MSH, XS, NPAR, BC, therm, hydro, xnp1, t+dt, dt, IV, resn, resnp1, epsi);
            if(NPAR.PreCond)
                [delta_vec]=OS_prec_tr(MSH, NPAR, XS, therm, hydro, BC, delta_vec, xnp1, t+dt, dt, IV);
            end
            iter_=iter_+iter;
            xnp1 = xnp1 + NPAR.damp*delta_vec;
            [resnp1,err2,therm] = comp_residual_tr( MSH, XS, NPAR, BC, therm, hydro, xnp1, t+dt, dt, IV, resn);
            % compute errors
            err=norm(delta_vec,inf);
            % printouts
            disp(sprintf('      JFNK iteration %i, gmres iter %i - %i,  error %g %g, tol %g',i,iter(1:2),err,err2,tol));
            %             if( err2 <1e-4 | i==NPAR.MaxNewton ), pause, end
            %             if( err2 <1e-9)
            if( err2 <10*tol)
                [Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xnp1 );
                [Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xnp1);
                subplot(2,1,1); plot(Teff);
                subplot(2,1,2); plot(Tmod);
                disp('converged nonlinear matrix-free');
                break
            else
                [Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xnp1 );
                [Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xnp1);
                subplot(2,1,1); plot(Teff);hold on;
                subplot(2,1,2); plot(Tmod);hold on;
                %                 if(i==1),pause,end
            end
        end
        Plevel(it+1)=therm.Plevel(end)
        disp('number of iterations');
        iter_
        disp('nbr time steps')
        nt
    else

        for os_iter=1:max_os_iter

            if(max_os_iter>1)
                aux = xnp1;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(neutro_first)
                if(do.neutro>0)
                    [Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xnp1 );
                    [Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xnp1);
                    [Anp1,ixs]=create_mat_tn(MSH,XS,NPAR,BC,Teff,Rmod,t+dt,IV);
                    % solve
                    xnp1(1:ns)=inv(( speye(ns)-dt/2*Anp1(1:ns,1:ns) )) ...
                        * ( speye(ns)+dt/2*An(1:ns,1:ns) ) * xn(1:ns);
                end
                [Pow,Plevel(it+1)] = comp_power( ixs, MSH, NPAR, xnp1(1:NPAR.neu), XS.G );
                if(test.thermo>0)
                    Plevel(it+1)=t_ramp(it+1);
                end
                [Pow,Plevel(it+1)] = norma_power( Pow, Plevel(it+1)*therm.Watt0, true );

                if(do.thermo>0)
                    [xnp1(NPAR.prec+1:NPAR.ther)] = solve_tr_hc( MSH, therm, hydro,...
                        xnp1(NPAR.prec+1:NPAR.ther),  xn(NPAR.prec+1:NPAR.ther), dt, xnp1(NPAR.ther+1:NPAR.siz), ...
                        xn(NPAR.ther+1:NPAR.siz), Pow, Pown, Tmod, NPAR.MatFree, NPAR.PreCond, 1, NPAR.method );

                end

                if(do.hydro)
                    if(test.hydro>0)
                        flxth=flxth0*t_ramp(it+1);
                    else
                        %%%%%%%% extract heat flux
                        [flxth]=extract_flxth( MSH, hydro, xnp1(NPAR.prec+1:NPAR.ther), Tmod, ...
                            xnp1(NPAR.ther+1:NPAR.siz) );
                    end
                    %%%%%%%% solve tr enthalpies
                    [xnp1(NPAR.ther+1:NPAR.siz)] = solve_tr_hy( MSH, therm, hydro, ...
                        xnp1(NPAR.ther+1:NPAR.siz), xn(NPAR.ther+1:NPAR.siz), dt, flxth, flxthn, ...
                        NPAR.MatFree, NPAR.PreCond, NPAR.method );

                    %%%%%%%% extract heat flux (to make sure consistency)
                    [flxth]=extract_flxth( MSH, hydro, xnp1(NPAR.prec+1:NPAR.ther), Tmod, ...
                        xnp1(NPAR.ther+1:NPAR.siz) );
                end
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(do.thermo>0)
                    [xnp1(NPAR.prec+1:NPAR.ther)] = solve_tr_hc( MSH, therm, hydro,...
                        xnp1(NPAR.prec+1:NPAR.ther),  xn(NPAR.prec+1:NPAR.ther), dt, xnp1(NPAR.ther+1:NPAR.siz), ...
                        xn(NPAR.ther+1:NPAR.siz), Pow, Pown, Tmod, NPAR.MatFree, NPAR.PreCond, 1, NPAR.method );

                end

                if(do.hydro)
                    if(test.hydro>0)
                        flxth=flxth0*t_ramp(it+1);
                    else
                        %%%%%%%% extract heat flux
                        [flxth]=extract_flxth( MSH, hydro, xnp1(NPAR.prec+1:NPAR.ther), Tmod, ...
                            xnp1(NPAR.ther+1:NPAR.siz) );
                    end
                    %%%%%%%% solve tr enthalpies
                    [xnp1(NPAR.ther+1:NPAR.siz)] = solve_tr_hy( MSH, therm, hydro, ...
                        xnp1(NPAR.ther+1:NPAR.siz), xn(NPAR.ther+1:NPAR.siz), dt, flxth, flxthn, ...
                        NPAR.MatFree, NPAR.PreCond, NPAR.method );

                    %%%%%%%% extract heat flux (to make sure consistency)
                    [flxth]=extract_flxth( MSH, hydro, xnp1(NPAR.prec+1:NPAR.ther), Tmod, ...
                        xnp1(NPAR.ther+1:NPAR.siz) );
                end

                if(do.neutro>0)
                    [Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xnp1 );
                    [Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xnp1);
                    [Anp1,ixs]=create_mat_tn(MSH,XS,NPAR,BC,Teff,Rmod,t+dt,IV);
                    % solve
                    xnp1(1:ns)=inv(( speye(ns)-dt/2*Anp1(1:ns,1:ns) )) ...
                        * ( speye(ns)+dt/2*An(1:ns,1:ns) ) * xn(1:ns);
                end
                [Pow,Plevel(it+1)] = comp_power( ixs, MSH, NPAR, xnp1(1:NPAR.neu), XS.G );
                [Pow,Plevel(it+1)] = norma_power( Pow, Plevel(it+1)*therm.Watt0, true );

            end

            % test for convergence
            if(max_os_iter>1)
                err = norm(aux-xnp1,inf);
                if(err<1E-8)
                    break
                elseif(os_iter==max_os_iter)
                    disp(sprintf('remaining error %g',err));
                    disp(sprintf('%i iterations',os_iter));
                    disp('not enough iterations, pb not CVed in os main')
                end
            end
        end

    end

    % output
    fwrite(fid,xnp1,'double');
    t=t+dt;

    xn     = xnp1;
    if(NPAR.MatFree==false)
        if(do.neutro>0)
            An     = Anp1;
        end
        Pown   = Pow;
        if(do.hydro>0)
            flxthn = flxth;
        end
        Tmodn  = Tmod;
    end
end

pp = Plevel/Plevel(1);
fwrite(fid,pp,'double');
fclose(fid)

figure(99);
tt=linspace(0,tend,nt+1);
figure(1);
plot(tt,pp,'r+-')
max(pp)

% % verify
% [Teff]       = extract_Teff( MSH, NPAR.prec, therm.rw, xn );
% [Rmod, Tmod] = extract_Rmod( MSH, NPAR.ther, hydro.pressure, xn);
% loc=[18,1]*1;
% [P,M,ixs] = create_mat_ss( MSH, XS, NPAR, BC, Teff, Rmod,loc,t);
% %%%%%%%%% solve ss neutronics
% [keff,vec(1:NPAR.neu),XS] = solve_ss(P,M,XS,NPAR.renorm);
% disp(sprintf('ss eigenvalue is %g',keff));
%
% rr=1-1/keff
% 649/(649-rr*1e5)
%
%
% %     save 1dresults_1sec_10ref.mat tt pp sav;
