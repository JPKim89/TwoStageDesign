[MM_D Opt_D] = TwoStageOptD(P0, P1, alfa, pwr)

%   This function computes sample size for two-stage design 
%   when the endpoint is response (binary) 

%   Inputs:
%   P0: the response rate under null hypothesis
%   P1: the response rate under alternative hypothesis
%   alfa and pwr: targeted type I error and statistical power

%   Outcomes:
%   MM_D: 	the results of Modified minimax design
%   Opt_D: 	the results of Modified Optimal design by Kim

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

If # of responders is less than equal to R1 in N1 subjects, then the study will be early terminated. 
If R1+1 or more responses are observed, then the second stage is open and (N-N1) additional subjects will be accrued.
If R+1 are noted in N subjects, the null hypothesis is rejected and the regimen seems promising for further study.
The expected sample size under the null hypothesis is EN. 
The probabilities of early stopping under null and alternative are PET0 and PET1, respectively.