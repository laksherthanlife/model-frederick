clear all

%% INITIALIZE TIME
to=0;
tf=120;
inc=120;
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

%% Phospho-tuning function: Change in Beta vs. time due to Hsf1 phosphorylation state (reflecting experimental data)
beta_lo = 1.5   ;
beta_hi = 10^4  ;
beta_IN = (beta_hi-beta_lo)./(1 + beta_lo*exp(-0.25*(time-60))) + beta_lo;   % WT phospho-tuning function 
beta_IN_dPO4 = ones(1,inc*1)*beta_lo;                                        % Delta PO4 = constant Beta

%% INITIAL CONDITIONS
HSPo	  = 1;          % Free HSP
Hsf1o     = 0;          % Free Hsf1
HSP_Hsf1o = 1/500;      % HSP-Hsf1 complex
HSP_UPo   = 0;          % HSP-UP complex
YFPo      = 3;          % YFP
UPo_39    = 10;
        
% Initial Condition Vector
IC = [HSPo, Hsf1o, UPo_39, HSP_Hsf1o, HSP_UPo, YFPo];

%% RUN ODEs
[t ,yWT] = ode23s(@(t,y)titration_YFP_FB_Phosph(t,y,k1,k2,k3,k4,k5,beta_IN,Kd,time), time, IC);  % 35°C ODE
[t ,ydPO4] = ode23s(@(t,y)titration_YFP_FB_Phosph(t,y,k1,k2,k3,k4,k5,beta_IN_dPO4,Kd,time), time, IC);  % 35°C ODE

%% PULL DATA
YFP_WT = yWT(:,6);      % WT Hsf1
YFP_dPO4 = ydPO4(:,6);  % Delta PO4 Hsf1 Mutant

%% PLOT

% Plot YFP vs Time for WT and Delta PO4 Mutant
figure
    plot(time,log2(YFP_WT/YFPo),'k'); hold on
    plot(time,log2(YFP_dPO4/YFPo),'r'); hold on
    set(gca,'FontSize',18)
    xlabel('Time (min)');
    ylabel('Foldchange YFP');
    xlim([0 time(end)])
    ylim([-.5 4])

% Plot Phospho-Tuning Function
figure
    semilogy(time,beta_IN,'k'); hold on
    semilogy(time,beta_IN_dPO4,'r'); hold on
    set(gca,'FontSize',18)
    xlabel('Time (min)');
    ylabel('Beta');
    xlim([0 time(end)])
    ylim([0.5 10^5])