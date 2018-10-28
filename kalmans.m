function [ xst_N Psk_N Lt] = kalmans(xpred,vpred,xfilt,vfilt,paramA)
T = length(xpred);
xst_N = cell(1,T);
xst_N(1,T) = xfilt(1,T);
Psk_N = cell(1,T);
Psk_N(1,T) = vfilt(1,T);
Lt = cell(1,T);
Lt(1,T) = {vfilt{1,T}*paramA'*inv(vpred{1,T})};

for i=T-1:-1:1
    Lt(1,i) = {vfilt{1,i}*paramA'*inv(vpred{1,i})};
    xst_N(1,i) = {xfilt{1,i} + Lt{1,i}*(xst_N{1,i+1}-xpred{1,i})};
    Psk_N(1,i) = {vfilt{1,i} + Lt{1,i}*(Psk_N{1,i+1}-vpred{1,i})*Lt{1,i}'};
end

end