clear all

%% INITIALIZE TIME
to=0;
tf=60;
inc=1001;
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


%% INITIAL CONDITIONS for Temp = 25°C
HSPo	  = 1;          % Free HSP
Hsf1o     = 0;          % Free Hsf1
UPo       = 0;          % Unfolded Protein
HSP_Hsf1o = 1/500;      % HSP-Hsf1 complex
HSP_UPo   = 0;          % HSP-UP complex
YFPo      = 0;          % YFP
ICo = [HSPo, Hsf1o, UPo, HSP_Hsf1o, HSP_UPo, YFPo];


%% INITIAL CONDITIONS for Temp = 35°C
HSPo	  = 1;          % Free HSP
Hsf1o     = 0;          % Free Hsf1
UPo       = 10;         % Unfolded Protein
HSP_Hsf1o = 1/500;      % HSP-Hsf1 complex
HSP_UPo   = 0;          % HSP-UP complex
YFPo      = 0;          % YFP
IC1 = [HSPo, Hsf1o, UPo, HSP_Hsf1o, HSP_UPo, YFPo];

%% RUN ODEs
[t ,y25] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, ICo);  % 25°C ODE
[t ,y35] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC1);  % 35°C ODE

%% PULL DATA
HSP_Hsf1_25 = y25(:,4);
HSP_Hsf1_35 = y35(:,4);
HSP_Hsf1_IPMS = HSP_Hsf1_35./HSP_Hsf1_25;      % Normalize Data

%% PLOT
figure
    plot(time,HSP_Hsf1_IPMS)
    set(gca,'FontSize',18)
    xlabel('Time (min)');
    ylabel('Normalized [HSP-Hsf1]');
    xlim([0 time(end)])
    ylim([0 1.2])
    