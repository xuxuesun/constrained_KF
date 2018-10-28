function [xt_t xtp_t Pt_t Ptp_t avgloglik] = kalmanf(sig,A,B,u,Q,R,C)
%state space model, e1~N(0,Q),e2~N(0,R)
%x = Ax+Bu+e1
%y = Cx+e2
xt0 = sig(:,1);
Pt0 = R;
T = size(sig,2);
xt_t = cell(1,T);
xt_t(1,1) = {xt0};
xtp_t = cell(1,T);
xtp_t(1,1) = {xt0};
resi = cell(1,T);
resi(1,1) = {[0
    0]};
Ptp_t = cell(1,T);
Ptp_t(1,1)={Pt0};
Pt_t = cell(1,T);
Pt_t(1,1) = {Pt0};
Rc = cell(1,T);
Rc(1,1) = {2*R};
K = cell(1,T);
K(1,1) = {[1
    1]};
loglik = [];

for i=2:length(sig)
    meas = sig(:,i);
    
    %prediction
    xtp_t(1,i) = {A*xt_t{1,i-1}+B*u};
    Ptp_t(1,i) = {A*Pt_t{1,i-1}*A'+Q};
    %residual
    resi(1,i) = {meas-C*xtp_t{1,i}};
    %residual covariance
    Rc(1,i) = {C*Ptp_t{1,i}*C'+R};
    %gain
    K(1,i) = {Ptp_t{1,i}*C'*inv(Rc{1,i})};
    %log likelihood   
    loglik = [loglik -T/2*log(2*pi)-log(det(Rc{1,i}))/2-resi{1,i}'*(Rc{1,i}\resi{1,i})/2];
    
    %update
    xt_t(1,i) = {xtp_t{1,i}+K{1,i}*resi{1,i}};
    Pt_t(1,i) = {Ptp_t{1,i}-K{1,i}*C*Ptp_t{1,i}};
end

avgloglik = sum(loglik)/length(loglik);
end