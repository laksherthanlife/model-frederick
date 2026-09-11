clear all

%% INITIALIZE TIME
to=0;
tf=120;
inc=121;
time = linspace(to,tf,inc);

%% ODE KINETIC/TXN PARAMETERS
% Obtained from Parameter Screen
k1 = 166.8;     % HSP-UP association
k2 =  2.783;    % HSP-Hsf1 disassociation
k3 = k1;        % HSP-Hsf1 association
k4 = 0.0464;    % HSP-Hsf1 disassociation
k5 = 4.642e-7;  % HSP-UP refolding into FP
beta = 1.7783;  % Txn activation
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
ICo = [HSPo, Hsf1o, UPo_25, HSP_Hsf1o, HSP_UPo, YFPo];
IC1 = [HSPo, Hsf1o, UPo_35, HSP_Hsf1o, HSP_UPo, YFPo];
IC2 = [HSPo, Hsf1o, UPo_39, HSP_Hsf1o, HSP_UPo, YFPo];
IC3 = [HSPo, Hsf1o, UPo_43, HSP_Hsf1o, HSP_UPo, YFPo];

%% RUN ODEs
[t ,y25] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, ICo);  % 25°C ODE
[t ,y35] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC1);  % 35°C ODE
[t ,y39] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC2);  % 39°C ODE
[t ,y43] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC3);  % 43°C ODE

%% PULL DATA
YFP_25 = y25(:,6);
YFP_35 = y35(:,6);
YFP_39 = y39(:,6);
YFP_43 = y43(:,6);

%% PLOT
   
figure
    plot(time,log2(YFP_35/YFPo),'k'); hold on
    plot(time,log2(YFP_39/YFPo),'r'); hold on
    plot(time,log2(YFP_43/YFPo),'b');
    set(gca,'FontSize',18)
    xlabel('Time (min)');
    ylabel('Foldchange YFP');
    xlim([0 time(end)])
    ylim([-.5 4])