clear all

%% EFFECT of lowering HSP70-Hsf1 affinity on heat shock time course

%% INITIALIZE TIME
to=0;
tf=240;
inc=241;
time = linspace(to,tf,inc);

%% ODE KINETIC/TXN PARAMETERS
% Obtained from Parameter Screen
k1 = 166.8;     % HSP-UP association

k2 =  2.783;    % HSP-Hsf1 disassociation
k2_low =  2.783*50;    % HSP-Hsf1 LOWER AFFINITY MUTANT

k3 = k1;        % HSP-Hsf1 association
k4 = 0.0464;    % HSP-Hsf1 disassociation
k5 = 4.642e-7;  % HSP-UP refolding into FP
beta = 1.7783/5;  % Txn activation (MODIFIED FROM FIRST PAPER)
Kd = 0.0022;    % Hsf1-DNA binding affinity


%% INITIAL CONDITIONS
HSPo	  = 1;          % Free HSP
Hsf1o     = 0;          % Free Hsf1
HSP_Hsf1o = 1/500;      % HSP-Hsf1 complex
HSP_UPo   = 0;          % HSP-UP complex
YFPo      = 3;          % YFP

% UP at different temperatures (From Temp vs UP relationship)
UPo_25 = 0.5183;
UPo_35 = 4.4492;
UPo_39 = 10.5141;
UPo_43 = 20.0397;         

% Initial Condition Vector
IC1 = [HSPo, Hsf1o, UPo_25, HSP_Hsf1o, HSP_UPo, YFPo];
IC2 = [HSPo, Hsf1o, UPo_39, HSP_Hsf1o, HSP_UPo];

%% RUN ODEs
    

% WT ODEs
[t ,y25] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC1);  % 25°C ODE
basal_WT = y25(end,6); % Get basal value of YFP reporter by taking steady-state value at 25°C

[t ,y39] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, [IC2 basal_WT]);  % 39°C ODE
YFP_39 = y39(:,6);

% Low affinity mutant ODE
[t ,y25] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2_low,k3,k4,k5,beta,Kd), time, IC1);  % 25°C ODE
basal_low = y25(end,6);

[t ,y39] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2_low,k3,k4,k5,beta,Kd), time, [IC2 basal_low]);  % 39°C ODE
YFP_39_low = y39(:,6);

%% PLOT
   
figure
    plot(time,YFP_39,'k'); hold on
    plot(time,YFP_39_low,'b'); hold on
    set(gca,'FontSize',18)
    xlabel('time (min)');
    ylabel('YFP (a.u.)');
    legend('WT','low affinity mutant')
    xticks([0,60,120,180,240])
    xlim([0 time(end)])
