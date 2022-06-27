clc;
clear;
RFD=@(x)ResidualFullData(x);
% initial approximation:
x0 = [1 1 1 1 1 1 1 1 1]; 
%'FinDiffType', 'central' вставить после 0.01
options = optimset('Display','Iter','LargeScale','on','MaxIter',20,...
         'TolX',1e-4,'DiffMinChange',1e-2, 'InitTrustRegionRadius', 0.01, 'FinDiffType', 'central', 'Algorithm', 'trust-region-reflective')                  %#ok<NOPRT>
lb = [0 0 0 0 0 0 0 0 0];
ub = [5 5 5 5 5 5 5 5 5];
% lb=[];
% ub=[];
% optimize
[x,resnor,resid,exit,output,lambda,jacob] = lsqnonlin(RFD,x0,lb,ub,options)
