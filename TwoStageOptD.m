function [MM_D Opt_D] = TwoStageOptD(P0, P1, alfa, pwr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function computes sample size for two-stage design 
%   when the endpoint is response (binary) 

%   Inputs:
%   P0: the response rate under null hypothesis
%   P1: the response rate under alternative hypothesis
%   alfa and pwr: targeted type I error and statistical power

%   Outcomes:
%   MM_D:  the results of Modified minimax design
%   Opt_D: the results of Modified Optimal design

%   For example,
%   MM_D.R1    : threshold for early termination
%   MM_D.N1    : Sample size of 1st stage
%   MM_D.R     : maximum number of responders to fail to reject p0
%   MM_D.N     : sample size of combined stage
%   MM_D.alpha : Actual alpha of minimax design
%   MM_D.power : Actual power of minimax design
%   MM_D.EN    : Expected sample size of combined stage under p0
%   MM_D.PET0  : Probability of Early Termination under p0
%   MM_D.PET1  : Probability of Early Termination under p1

%   This program calls either "Modified Two-Stage Design database.mat" or
%   "TwoStageDesign.m"

% If # of responders is less than equal to R1 in N1 subjects, then the study will be early terminated. 
% If R1+1 or more responses are observed, then the second stage is open and (N-N1) additional subjects will be accrued.
% If R+1 are noted in N subjects, the null hypothesis is rejected and the regimen seems promising for further study.
% The expected sample size under the null hypothesis is EN. 
% The probabilities of early stopping under null and alternative are PET0 and PET1, respectively.

%   Coded by Jongphil Kim on April 14, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   1. Check if the solution is available from the database
load('Modified Two-Stage Design database.mat');
ID_MM = (MiniMaxD(:,1) == alfa) & (MiniMaxD(:,2) == pwr) & ...
    (MiniMaxD(:,3) == P0) & (MiniMaxD(:,4) == P1 ) ;
MM_sol = MiniMaxD(ID_MM, 5:end);

if isempty(MM_sol)
    ID_MM = (abs(MiniMaxD(:,1) - alfa) < eps) & (abs(MiniMaxD(:,2) - pwr) < eps) & ...
    (abs(MiniMaxD(:,3) - P0) < eps) & (abs(MiniMaxD(:,4) - P1) < eps) ;
    MM_sol = MiniMaxD(ID_MM, 5:end);

    ID_Opt = (abs(OptimalD(:,1) - alfa) < eps) & (abs(OptimalD(:,2) - pwr) < eps) & ...
        (abs(OptimalD(:,3) - P0) < eps) & (abs(OptimalD(:,4) - P1) < eps) ;
    Opt_sol = OptimalD(ID_Opt, 5:end);
else
    ID_Opt = (OptimalD(:,1) == alfa) & (OptimalD(:,2) == pwr) & ...
        (OptimalD(:,3) == P0) & (OptimalD(:,4) == P1) ;
    Opt_sol = OptimalD(ID_Opt, 5:end);
end

%   2. Find the optimal solution
if ~isempty(MM_sol)    
    MM_D.R1 = MM_sol(1,1);          MM_D.N1 = MM_sol(1,2);
    MM_D.R  = MM_sol(1,3);          MM_D.N  = MM_sol(1,4);
    MM_D.alfa = MM_sol(1,5);        MM_D.power = MM_sol(1,6);
    MM_D.EN = MM_sol(1,7);          MM_D.PET0 = MM_sol(1,8);
    MM_D.PET1 = MM_sol(1,9);

    Opt_D.R1 = Opt_sol(1,1);        Opt_D.N1 = Opt_sol(1,2);
    Opt_D.R  = Opt_sol(1,3);        Opt_D.N  = Opt_sol(1,4);
    Opt_D.alfa = Opt_sol(1,5);      Opt_D.power = Opt_sol(1,6);
    Opt_D.EN = Opt_sol(1,7);        Opt_D.PET0 = Opt_sol(1,8);
    Opt_D.PET1 = Opt_sol(1,9);
else
    [MM_D Opt_D] = TwoStageDesign(P0, P1, alfa, pwr, 1);
end


disp('************************************************************************************');
disp(['P0 = ', num2str(P0, '%3.2f'), ', P1 = ', num2str(P1, '%3.2f'), ...
    ', alpha = ', num2str(alfa, '%3.2f'), ', and power = ', num2str(pwr, '%3.2f'), '.']);  

disp('************************************************************************************');
disp('Modified Minimax Two-Stage Design Results:');
disp(['r1 = ', num2str(MM_D.R1, '%3.0f'), ', n1 = ', num2str(MM_D.N1, '%3.0f'), ...
         ', r = ', num2str(MM_D.R, '%3.0f'), ', and n = ', num2str(MM_D.N, '%3.0f'), '.']);    
disp(['Type I error = ', num2str(MM_D.alfa, '%5.4f'), ', Statistical Power = ',...
        num2str(MM_D.power, '%5.4f')]);
disp('If the null hypothesis is true, then');
disp(['the expected sample size = ', num2str(MM_D.EN, '%5.2f'), ...
        ', the probability of early termination = ', num2str(MM_D.PET0, '%5.3f'), ' under the null.']);
disp(['and the probability of early termination = ', num2str(MM_D.PET1, '%5.3f'), ' under the alternative.']);
disp('************************************************************************************');
disp('Modified Optimal Two-Stage Design Results:');
disp(['r1 = ', num2str(Opt_D.R1, '%3.0f'), ', n1 = ', num2str(Opt_D.N1, '%3.0f'), ...
         ', r = ', num2str(Opt_D.R, '%3.0f'), ', and n = ', num2str(Opt_D.N, '%3.0f'), '.']);    
disp(['Type I error = ', num2str(Opt_D.alfa, '%5.4f'), ', Statistical Power = ',...
        num2str(Opt_D.power, '%5.4f')]);
disp('If the null hypothesis is true, then');
disp(['the expected sample size = ', num2str(Opt_D.EN, '%5.2f'), ...
        ', the probability of early termination = ', num2str(Opt_D.PET0, '%5.3f'), ' under the null.']);
disp(['and the probability of early termination = ', num2str(Opt_D.PET1, '%5.3f'), ' under the alternative.']);
disp('************************************************************************************');    

end

