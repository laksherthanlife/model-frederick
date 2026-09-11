clear all

%% INITIALIZE TIME
to=0;
tf=120;
inc=121;
time = linspace(to,tf,inc);

%% HSF1 TITRATION RANGE
Hsf1_add = logspace(-4,3,50); 

%% ODE KINETIC/TXN PARAMETERS
% Obtained from Parameter Screen
k1 = 166.8;     % HSP-UP association
k2 =  2.783;    % HSP-Hsf1 disassociation
k3 = k1;        % HSP-Hsf1 association
k4 = 0.0464;    % HSP-Hsf1 disassociation
k5 = 4.642e-7;  % HSP-UP refolding into FP
beta = 1.7783;  % Txn activation
Kd = 0.0022;    % Hsf1-DNA binding affinity

for i=1:length(Hsf1_add)

    %% Initial Conditions for Est-Hsf1(WT)
    HSPo	  = 1;              % Free HSP
    Hsf1o     = 0;              % Free Hsf1
    UPo       = 0;              % Unfolded Protein
    HSP_Hsf1o = Hsf1_add(i);    % HSP-Hsf1 complex
    HSP_UPo   = 0;              % HSP-UP complex
    YFPo      = 3;              % Initial YFP concentration
    ICo = [HSPo, Hsf1o, UPo, HSP_Hsf1o, HSP_UPo, YFPo];
    
    %% Initial Conditions for Est-Hsf1(?DBD)
    HSPo	  = 1;              % Free HSP
    Hsf1o     = 0;              % Free Hsf1
    UPo       = 0;              % Unfolded Protein
    HSP_Hsf1o = 1/500;          % HSP-Hsf1 complex
    HSP_UPo   = 0;              % HSP-UP complex
    YFPo      = 3;              % Initial YFP concentration
    Hsf1_decoy = Hsf1_add(i);	% Hsf1 Decoy concentration
    HSP_Hsf1_decoy = 0;         % HSP_Hsf1_decoy concentraiton
    
    IC1 = [HSPo, Hsf1o, UPo, HSP_Hsf1o, HSP_UPo, YFPo, Hsf1_decoy, HSP_Hsf1_decoy];
    
    %% Run the ODEs
    [t ,y] = ode23s(@(t,y)titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd), time, ICo);
    [t2 ,y2] = ode23s(@(t,y)titration_YFP_FB_decoy(t,y,k1,k2,k3,k4,k5,beta,Kd), time, IC1);

    % WT Hsf1 Expression
    Hsf1_free(i) = y(end,2);         % Steady state of free Hsf1
    YFP(i) = y(end,6);               % Steady state of YFP
    
    % Decoy Hsf1 Expression
    Hsf1_free2(i) = y2(end,2);       % Steady state of free Hsf1
    YFP2(i) = y2(end,6);             % Steady state of YFP
end

%% PLOT

figure
    semilogx(Hsf1_add,YFP,'k'); hold on;
    semilogx(Hsf1_add,YFP2,'r');
    set(gca,'FontSize',18)
    xlabel('[HSF1] Titration');
    ylabel('YFP (a.u.)');
    xlim([Hsf1_add(1) Hsf1_add(end)])   
    
figure
    semilogx(Hsf1_add,log2(YFP/YFPo),'k'); hold on;
    semilogx(Hsf1_add,log2(YFP2/YFPo),'r');
    set(gca,'FontSize',18)
    xlabel('[HSF1] Titration');
    ylabel('Foldchange YFP');
    xlim([Hsf1_add(1) Hsf1_add(end)])   
    ylim([-1 5])