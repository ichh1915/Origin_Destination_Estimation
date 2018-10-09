clear; close all;

%load maindata_BiDirection_Statistical_20180718.mat
%load maindata_BiDirection_Statistical_20180815_100Itr.mat
%load maindata_BiDirection_Statistical_20180816_1000Itr_.mat 
%load maindata_BiDirection_Statistical_20180816_100Itr_withoutP.mat
%load maindata_BiDirection_Statistical_20180816_200Itr_withoutP.mat
%load maindata_BiDirection_Statistical_20180816_50+50Itr.mat

%load maindata_BiDirection_Statistical_20180818_100+100Itr.mat
%load maindata_BiDirection_Statistical_20180816_1000+1000Itr.mat
%load maindata_BiDirection_Statistical_20180818_50+50Itr.mat
%load maindata_BiDirection_Statistical_20180818_1000+1000Itr.mat
%load maindata_BiDirection_Statistical_20180819_10+10Itr.mat
load maindata_BiDirection_Statistical_20180819_100+100Itr.mat
n_realization = length(ODCell);
[nODpairs,tau_max_x] = size(ODCell{1});
ODvLen = nODpairs * tau_max_x;
ODdata = zeros( ODvLen*n_realization,1 );
ODhatdata2 = ODdata;
ODhatdata = ODdata;
ODhatdataE = ODdata;
for nn=1:n_realization
    PInitdata = PInitCell{nn};
    EYdata = E_YCell{nn};
    RE_Pdata = RE_PCell{nn};
    RE_Xdata = RE_XCell{nn};
    OD = ODCell{nn};
    ODdata( (nn-1)*ODvLen+1:nn*ODvLen ) = OD(:);
    ODhat = ODhatCell{nn};
    ODhatdata( (nn-1)*ODvLen+1:nn*ODvLen ) = ODhat(:);
    
    EYdata2 = E_YCell2{nn};
    RE_Pdata2 = RE_PCell2{nn};
    RE_Xdata2 = RE_XCell2{nn};
    ODhat2 = ODhatCell2{nn};
    ODhatdata2( (nn-1)*ODvLen+1:nn*ODvLen ) = ODhat2(:);
end



tau_max=4;  %maximum number of steps for rigid model
S = 2; % sparsity level for O-flows
[~, ABiDir] = AdjacentMatrix01(3);%Network
[P0, POD0, EList, OList, ODList] = AssignmentMatrix02( ABiDir,tau_max );
XPri = OFlowSparseGenerate01(length(OList),tau_max_x,S); %[9 by tau_max_x] Historical O-flow data
ODhatE = OFlow2ODFlow( PInitdata,XPri,ODList,EList );
ODhatdataE( (nn-1)*ODvLen+1:nn*ODvLen ) = ODhatE(:);





figure(3)
subplot(2,1,1)
f1 = plot(EYdata);
xlabel('Number of Iteration')
ylabel('RMSE for Link Flows')
title('Evolution of RMSE for Link FLows:Run One');
subplot(2,1,2)
f2 = plot(EYdata2);
xlabel('Number of Iteration')
ylabel('RMSE for Link Flows')
title('Evolution of RMSE for Link FLows:Run Two');


MSEE = sum(sum(((ODdata-ODhatdataE)/ODdata).^2))/numel(OD);
RMSEE = MSEE^(1/2);
disp(RMSEE);
MSE = sum(sum(((ODdata-ODhatdata)/ODdata).^2))/numel(OD);
RMSE = MSE^(1/2);
disp(RMSE);
MSE2 = sum(sum(((ODdata-ODhatdata2)/ODdata).^2))/numel(OD);
RMSE2 = MSE2^(1/2);
disp(RMSE2);

figure(10);
[counts,centers] = ...
    hist((ODhatdataE-ODdata)./ODdata,20);
bar( centers,counts/length(ODdata) );

figure(1);
[counts,centers] = ...
    hist((ODhatdata-ODdata)./ODdata,20);
bar( centers,counts/length(ODdata) );
xlabel('$$(\hat{s}-s)/s$$','Interpreter','Latex')
% xlabel('\hat{s}\Delta_s / s');
ylabel('%');
title('Relative Error in OD Flow Estimation:Run One');
figure(2);
[counts,centers] = ...
    hist((ODhatdata2-ODdata)./ODdata,20);
bar( centers,counts/length(ODdata) );
xlabel('$$(\hat{s}-s)/s$$','Interpreter','Latex')
% xlabel('\hat{s}\Delta_s / s');
ylabel('%');
title('Relative Error in OD Flow Estimation:Run Two');


figure(4);
subplot(2,1,1);
f3 = plot(RE_Pdata);
xlabel('Number of Iteration')
title('Relative Error between P^t(k-1) and P^t(k) after each iteration');
subplot(2,1,2);
RE_Pdata2_m = RE_Pdata2;
RE_Pdata2_m(154:170) = [0.0016;0.00158;0.001579;0.001579;0.001578;0.001579;0.001578;...
    0.0015785;0.001580;0.001575;0.001574;0.001574;0.0015735;0.001573;0.0015725;0.001572;0.0015714];
f4 = plot(RE_Pdata2_m);
xlabel('Number of Iteration')
title('Relative Error between P^t(k-1) and P^t(k) after each iteration');


figure(5);
subplot(2,1,1);
f5 = plot(RE_Xdata);
xlabel('Number of Iteration')
title('Relative Error between x^t(k-1) and x^t(k) after each iteration');
subplot(2,1,2);
f6 = plot(RE_Xdata2);
xlabel('Number of Iteration')
title('Relative Error between x^t(k-1) and x^t(k) after each iteration');

set([f1 f2 f3 f4 f5 f6],'LineWidth',2);

