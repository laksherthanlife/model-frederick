function ty=titration_YFP_FB(t,y,k1,k2,k3,k4,k5,beta,Kd)

%% DEFINE PARAMETERS

% Txn Kinetics
nH = 3;     % Hsf1 Trimer
kdil = 0;   % Approximate to 0 because of short experiments

%% DEFINE VARIABLES

HSP	= y(1);         % Free HSP
Hsf1 = y(2);        % Free Hsf1
UP = y(3);          % Unfolded Protein
HSP_Hsf1 = y(4);    % HSP-Hsf1 complex
HSP_UP = y(5);      % HSP-UP complex
YFP = y(6);         % YFP


%% DIFFERENTIAL EQUATIONS
% Equations are based on simple mass action kinetics 
% The exceptions are transcriptional induction of HSPs and Reporter, which use Hill function

% d[HSP]/dt
ty(1) = k2*HSP_Hsf1 - k1*HSP*Hsf1 + k4*HSP_UP - k3*HSP*UP + k5*HSP_UP + beta*Hsf1^nH/(Kd^nH + Hsf1^nH);

% d[Hsf1]/dt
ty(2) = k2*HSP_Hsf1 - k1*HSP*Hsf1 ;

% d[UP]/dt
ty(3) = k4*HSP_UP - k3*HSP*UP - k5*HSP_UP;

% d[HSP-Hsf1]/dt
ty(4) = k1*HSP*Hsf1 - k2*HSP_Hsf1 ;

% d[HSP-UP]/dt
ty(5) = k3*HSP*UP - k4*HSP_UP - k5*HSP_UP;

% d[YFP]/dt
ty(6) = beta*Hsf1^nH/(Kd^nH + Hsf1^nH) - kdil*YFP;

% matrix equation to be evaluated by ODE solver 
ty = [ty(1);ty(2);ty(3);ty(4);ty(5);ty(6)];

end

 