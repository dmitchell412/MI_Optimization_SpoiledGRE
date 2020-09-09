
% MI-based optimization of parameter space using fminsearch
% MI calculated by Gauss-Hermite quadrature

function [MIobjfun]=MI_objfun_ernst_nonlin(flipAngle,tisinput,acqparam,signu)

%% Assign Acquisition Parameters
% Default Parameters
K=acqparam(1);
H=acqparam(2);
TR=acqparam(3);
TE=acqparam(4);

%% Generate Quadrature Points for MI Calculation
NumQP=5;
[~,xn,~,~,wn]=GaussHermiteNDGauss(NumQP,tisinput(1),tisinput(3));

%% Evaluate Signal Model
lqp=length(xn{1}(:));
parfor qp=1:lqp
    
    S(qp)=spoiledgre(K,H,flipAngle,TR,TE,xn{1}(qp),tisinput(2));
    
end

%% Gauss-Hermite Quadrature MI Approximation
[~,xn2,~,~,wn2]=GaussHermiteNDGauss(NumQP,0,signu);
S2=repmat(S,[size(xn2{1},1),1])+repmat(xn2{1},[1,size(S,2)]);

lnterm=log(sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2));
pterm=sum(repmat(wn(:)',[size(xn2{1},1),1]).*S2,2);
hz=-sum(wn2.*pterm.*lnterm);
MI=hz;
MIobjfun=-real(MI);

end

function S = spoiledgre(K,H,flipAngle,TR,TE,T1,T2star)
    S = (K.*H.*sind(flipAngle).*(1-exp(-TR./T1)).*exp(-TE./T2star))./(1-exp(-TR./T1).*cosd(flipAngle));
end