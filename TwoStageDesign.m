function [MinimaxD OptimalD] = TwoStageDesign(pp0, pp1, T_alfa, T_PWR, NoShow, GammaValue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function computes sample size for two-stage design 
%   when the endpoint is response (binary) 

%   Inputs:
%   p0: the response rate under null hypothesis
%   p1: the response rate under alternative hypothesis
%   alfa and beta: targeted type I error and statistical power
%   if NoShow = 0               Display the results on screen by default
%   if NoShow = Other values    Do Not Display the results on screen
%   GammaValue = [gamma_1; gamma_2], gamma_1 = 1/3 and gamma_2 = 2/3

%   Outcomes:
%   MinimaxD: the results of Modified minimax design
%   OptimalD: the results of Modified Optimal design

%   For example,
%   MinimaxD.R1    : threshold for early termination
%   MinimaxD.N1    : Sample size of 1st stage
%   MinimaxD.R     : maximum number of responders to fail to reject p0
%   MinimaxD.N     : sample size of combined stage
%   MinimaxD.alpha : Actual alpha of minimax design
%   MinimaxD.power : Actual power of minimax design
%   MinimaxD.EN    : Expected sample size of combined stage under p0
%   MinimaxD.PET0  : Probability of Early Termination under p0
%   MinimaxD.PET1  : Probability of Early Termination under p1

%   Coded by Jongphil Kim, last update on May 08, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    error('4 arguments should be provided');
elseif nargin == 4
    NoShow = 0;             %   Display the results on screen by default
    Gamma_1 = 1/3;          %   Default GammaValue
    Gamma_2 = 2/3;
elseif nargin == 5
    Gamma_1 = 1/3;          %   Default GammaValue
    Gamma_2 = 2/3;
elseif nargin == 6
    GammaValue = GammaValue(:);
    Gamma_1 = GammaValue(1,1);
    Gamma_2 = GammaValue(2,1);       
elseif nargin > 6
    error('4 to 6 arguments should be provided');
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       1. Parameter Settings
NNmax = 150;                %   Arbitrarily chosen number for upper limit
Epsilon = 0.1;              %   PET_1 <= 0.1

StartN = floor(0.7*BinoSample (pp0, pp1, T_alfa, T_PWR, 1));
Ns = (StartN : NNmax)';     %   feasible set for maximum sample size
SearchMax = 15;             %   just for optimal design




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Minimax solution

disp('Find the modified minimax two-stage design');
Feasible_MM = [];
for ii = 1 : length(Ns)     %   For a given n
    
    NN = Ns(ii,1);
    if Gamma_1*NN > floor(Gamma_1*NN)                           % non-integer
        NN1 = (floor(Gamma_1*NN)+1 : floor(Gamma_2*NN))';       % constraint on n1
    else
        NN1 = (floor(Gamma_1*NN) : floor(Gamma_2*NN))';         % constraint on n1
    end    
    
    for jj = 1 : length(NN1)
        RR1 = (0 : NN1(jj, 1) - 1)';
        RR1 = RR1(binocdf(RR1, NN1(jj,1), pp1) <= Epsilon, 1); % PET1 <= epsilon        
        NN2 = NN - NN1(jj,1);
        
        if ~isempty(RR1)        
            for kk = 1 : length(RR1)
                PET_p0 = binocdf(RR1(kk,1), NN1(jj,1), pp0);
                PET_p1 = binocdf(RR1(kk,1), NN1(jj,1), pp1);
                
                RR = (RR1(kk,1) : RR1(kk,1) + NN2)';
                TP1 = TwoStageAlpha(RR1(kk,1), NN1(jj,1), RR, NN, pp0);
                Pwr1 = TwoStagePower(RR1(kk,1), NN1(jj,1), RR, NN, pp1);

                IDs = (TP1 <= T_alfa) & (Pwr1 >= T_PWR);    % constraints on alpha and beta
                if sum(IDs) >= 1
                    RR = RR(IDs , 1);   LRR = length(RR);
                    TP1 = TP1(IDs, 1);  Pwr1 = Pwr1(IDs, 1);
                    ExpN = NN1(jj, 1) + (1 - PET_p0)*NN2;
                    Feasible_MM = [Feasible_MM; repmat([RR1(kk,1) NN1(jj,1)], ones(LRR),1) ...
                        RR repmat(NN, LRR,1) TP1 Pwr1 repmat([ExpN PET_p0 PET_p1], ones(LRR),1)];
                else
                end
            end
        else           
        %   Move the next NN1            
        end
    end
    
    if ~isempty(Feasible_MM)
        break;
    else
    end
    
    disp([Ns(ii,1)]);
end

MinimaxD.feasible = Feasible_MM;    %   All solutions statisfying two constraints
Sol = sortrows(Feasible_MM, 7);     % find solution which minimizes EN

MinimaxD.R1   = Sol(1,1);       MinimaxD.N1   = Sol(1,2);
MinimaxD.R    = Sol(1,3);       MinimaxD.N    = Sol(1,4);
MinimaxD.alfa = Sol(1,5);       MinimaxD.power = Sol(1,6);
MinimaxD.EN   = Sol(1,7);       MinimaxD.PET0 = Sol(1,8);
MinimaxD.PET1 = Sol(1,9);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Optimal solution
disp('Find the modified Optimal two-stage design');

OptR1 = MinimaxD.R1;            OptN1 = MinimaxD.N1;
OptR = MinimaxD.R;              OptN =  MinimaxD.N;
OptAlfa = MinimaxD.alfa;        OptPower = MinimaxD.power;
OptEN = MinimaxD.EN;            OptPET0 = MinimaxD.PET0;
OptPET1 = MinimaxD.PET1;

FeasibleSet = [OptR1 OptN1 OptR OptN OptAlfa OptPower OptEN OptPET0 OptPET1];
Nmax = min(OptN + SearchMax, NNmax);      %   Upper limit for N in optimal design
N = OptN + 1;

while N <= Nmax 

    if Gamma_1*N > floor(Gamma_1*N)                         % non-integer
        N1 = (floor(Gamma_1*N)+1 : floor(Gamma_2*N))';      % constraint on n1
    else
        N1 = (floor(Gamma_1*N) : floor(Gamma_2*N))';        % constraint on n1
    end    
    disp([N Nmax]);
    
    for jj = 1 : length(N1)
        RR1 = (0 : N1(jj, 1) - 1)';
        N2 = N - N1(jj,1);
                
        for kk = 1 : length(RR1)
            PET_p0 = binocdf(RR1(kk,1), N1(jj,1), pp0);
            PET_p1 = binocdf(RR1(kk,1), N1(jj,1), pp1);
            ExpN = N1(jj, 1) + (1 - PET_p0)*N2;
            
            if (ExpN < OptEN) && (PET_p1 <= Epsilon)    % Optimal sol and constraints on PET1
                RR = (RR1(kk,1) : RR1(kk,1) + N2)'; 
            
                for ll = 1 : length(RR)
                    TP1 = TwoStageAlpha(RR1(kk,1), N1(jj,1), RR(ll,1), N, pp0);
                    
                    if TP1 <= T_alfa                    % Constraints on alpha
                        Pwr1 = TwoStagePower(RR1(kk,1), N1(jj,1), RR(ll,1), N, pp1);
                        
                        if Pwr1 >= T_PWR                % Constraints on power
                            FeasibleSet = [FeasibleSet; RR1(kk,1) N1(jj,1) RR(ll,1) N TP1 Pwr1 ExpN PET_p0 PET_p1];
                            Nmax = min(N+SearchMax, NNmax);     % Reset on the maximum of n
                            OptEN = ExpN;               % update on optimal solution
                        else
                            % Void
                        end
                    else
                        % Void
                    end
                end
            else
                % Void
            end
        end
    end
    N = N + 1;
end

OptimalD.Feasible = FeasibleSet;    Sol = sortrows(FeasibleSet, 7);
OptimalD.R1 = Sol(1,1);             OptimalD.N1 = Sol(1,2);
OptimalD.R  = Sol(1,3);             OptimalD.N = Sol(1,4);
OptimalD.alfa = Sol(1,5);           OptimalD.power = Sol(1,6);
OptimalD.EN = Sol(1,7);             OptimalD.PET0 = Sol(1,8);
OptimalD.PET1 = Sol(1,9);

if NoShow == 0
    disp('************************************************************************************');
    disp(['Gamma 1 (lower limit) = ', num2str(Gamma_1, '%4.3f'), ', Gamma 2 (upper limit) = ', num2str(Gamma_2, '%4.3f')]);
    disp('************************************************************************************');
    disp('Modified Minimax Two-Stage Design Results:');
    disp(['r1 = ', num2str(MinimaxD.R1, '%3.0f'), ', n1 = ', num2str(MinimaxD.N1, '%3.0f'), ...
         ', r = ', num2str(MinimaxD.R, '%3.0f'), ', and n = ', num2str(MinimaxD.N, '%3.0f'), '.']);    
    disp(['Type I error = ', num2str(MinimaxD.alfa, '%5.4f'), ', Statistical Power = ',...
        num2str(MinimaxD.power, '%5.4f')]);
    disp('If the null hypothesis is true, then');
    disp(['the expected sample size = ', num2str(MinimaxD.EN, '%5.2f'), ...
        ', the probability of early termination = ', num2str(MinimaxD.PET0, '%5.3f'), ' under the null.']);
    disp(['and the probability of early termination = ', num2str(MinimaxD.PET1, '%5.3f'), ' under the alternative.']);
    disp('************************************************************************************');
    disp('Modified Optimal Two-Stage Design Results:');
    disp(['r1 = ', num2str(OptimalD.R1, '%3.0f'), ', n1 = ', num2str(OptimalD.N1, '%3.0f'), ...
         ', r = ', num2str(OptimalD.R, '%3.0f'), ', and n = ', num2str(OptimalD.N, '%3.0f'), '.']);    
    disp(['Type I error = ', num2str(OptimalD.alfa, '%5.4f'), ', Statistical Power = ',...
        num2str(OptimalD.power, '%5.4f')]);
    disp('If the null hypothesis is true, then');
    disp(['the expected sample size = ', num2str(OptimalD.EN, '%5.2f'), ...
        ', the probability of early termination = ', num2str(OptimalD.PET0, '%5.3f'), ' under the null.']);
    disp(['and the probability of early termination = ', num2str(OptimalD.PET1, '%5.3f'), ' under the alternative.']);
    disp('************************************************************************************');    
else
end

end
% end of function




