clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t2,p2]=post_power_level_adapt('y_case_odetb_e2.mat');
[t3,p3]=post_power_level_adapt('y_case_odetb_e3.mat');
[t4,p4]=post_power_level_adapt('y_case_odetb_e4.mat');

[tOS_notcv_1,pOS_notcv_1]=post_power_level_fixed('Case1_OS_1it_5sec_cv_1.txt');
[tOS_cv_1,pOS_cv_1]=post_power_level_fixed('Case1_OS_100it_5sec_cv_1.txt');

[tOS_notcv_01,pOS_notcv_01]=post_power_level_fixed('Case1_OS_1it_5sec_cv_01.txt');
[tOS_cv_01,pOS_cv_01]=post_power_level_fixed('Case1_OS_100it_5sec_cv_01.txt');

figure(49);
plot(t2,p2,'r+-','LineWidth',2); hold on
plot(t3,p3,'ko-','LineWidth',2); hold on
plot(t4,p4,'bs-'); hold on
xlabel('time, sec','FontSize',12)
ylabel('Normalized Power','FontSize',12)
axis([0. 5 1 1.85]);
grid on;
legend('tol=1e-2 ','tol=1e-3 ','tol=1e-4 ');

r= max(p4)/max(pOS_cv_01)*1.007

figure(50);
plot(t3,p3,'bs-','LineWidth',2); hold on
plot(tOS_notcv_1,r*pOS_notcv_1,'rv-','LineWidth',2); hold on
plot(tOS_cv_1,r*pOS_cv_1,'kv-','LineWidth',2); hold on
xlabel('time, sec','FontSize',12)
ylabel('Normalized Power','FontSize',12)
grid on;
axis([0.1 1.5 1.2 1.95]);
legend('Time Adapt RK, tol=1e-3 ',...
       'Crank-Nicholson, no iter',...
       'Crank-Nicholson, FPI    ');

   
r= max(p4)/max(pOS_cv_01)*1.001

figure(51);
plot(t3,p3,'bs-','LineWidth',2); hold on
plot(tOS_notcv_01,r*pOS_notcv_01,'rv-','LineWidth',2); hold on
plot(tOS_cv_01,r*pOS_cv_01,'kv-','LineWidth',2); hold on
xlabel('time, sec','FontSize',12)
ylabel('Normalized Power','FontSize',12)
grid on;
axis([0.2 .5 1.6 1.85]);
legend('Time Adapt RK, tol=1e-3 ',...
       'Crank-Nicholson, no iter',...
       'Crank-Nicholson, FPI    ');
   
[max(p2) max(p3) max(p4); length(t2) length(t3) length(t4)]

[max(pOS_notcv_1) max(pOS_cv_1) max(pOS_notcv_01) max(pOS_cv_01);...
    length(tOS_notcv_1) length(tOS_cv_1) length(tOS_notcv_01) length(tOS_cv_01) ]
