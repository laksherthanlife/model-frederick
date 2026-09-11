clear all

%% TITRATING HSP70-HSF1 binding affinity AFFECTS YFP reporter basal activity

%% INITIALIZE TIME
to=0;
tf=240;
inc=241;
time = linspace(to,tf,inc);

%% ODE KINETIC/TXN PARAMETERS
% Newer parameters
k1 = 166.8;     % HSP-UP association
k2_r = logspace(0,4,20);    % HSP-Hsf1 disassociation range
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
IC2 = [HSPo, Hsf1o, UPo_39, HSP_Hsf1o, HSP_UPo, YFPo];

%% RUN ODEs
output = zeros(size(k2_r));

for i=1:length(k2_r)
    
    k2 = k2_r(i);   % HSP-Hsf1 disassociation
    [t ,y25] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC1);  % 25°C ODE
    output(i) = y25(end,6); % Take steady-state level of YFP reporter at 25°C
      
end

%% PLOT
   
figure
    semilogx(k2_r,output/YFPo)
    xlabel('Hsf1-HSP disassociation rate');
    ylabel('fold change in basal YFP');
    set(gca,'FontSize',18)
    xlim([min(k2_r) max(k2_r)])
    ylim([0 12])
