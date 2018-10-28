clc;
clear all;
iter_cnt = 0;
difQ = 100;
difR = 100;

load 'com17.mat';
kalfilt.sig = com(2,:);
kalfilt.a = 0;
kalfilt.b = 0;
Ny = size(kalfilt.sig,1);
Nx = Ny;
kalfilt.F =  eye(Nx);   %state transition model
kalfilt.H = eye(Ny);    %observation model
kalfilt.Qt = 0;
kalfilt.Rt = cov(kalfilt.sig(1,1:50));
convergen_flag = 0;
loglikval = 100;
while ~convergen_flag           %convergence (F,H not update)   
    %filter data
    
    loglikval0 = loglikval;
    Rt0 = kalfilt.Rt;
    Qt0 = kalfilt.Qt;
    
    [xtt xtpt Ptt Ptpt loglikval] = kalmanf(kalfilt.sig,kalfilt.F,0,0,kalfilt.Qt,kalfilt.Rt,kalfilt.H);

    %smooth data
    [xsm psm ltval] = kalmans(xtpt,Ptpt,xtt,Ptt,kalfilt.F);
        
    %update parameter
    tempval = 0;
    for i=2:length(xsm)
       tempval = tempval+(xsm{1,i}-kalfilt.F*xsm{1,i-1})*(xsm{1,i}-kalfilt.F*xsm{1,i-1})' + kalfilt.F*psm{1,i-1}*kalfilt.F'+psm{1,i}-psm{1,i}*ltval{1,i-1}'*kalfilt.F'-kalfilt.F*ltval{1,i-1}*psm{1,i};
    end
    kalfilt.Qt = tempval/length(xsm);
    tempval = 0;
    for i=1:size(kalfilt.sig,2)
        tempval = tempval+ (kalfilt.sig(i)-kalfilt.H*xsm{1,i})*(kalfilt.sig(i)-kalfilt.H*xsm{1,i})'+kalfilt.H*psm{1,i}*kalfilt.H';
    end
    kalfilt.Rt = tempval/size(kalfilt.sig,2);
    if norm(Qt0-kalfilt.Qt,Inf)<1e-6 && norm(Rt0-kalfilt.Rt,Inf)<1e-6 && abs(loglikval0-loglikval)<1e-6
        convergen_flag = 1;
    end
    
    [xtt2 xtpt2 Ptt2 Ptpt2 loglikval2] = kalmanf(kalfilt.sig,kalfilt.F,0,0,kalfilt.Qt,kalfilt.Rt,kalfilt.H);
    plot(kalfilt.sig(1,:),'.');hold on;
    for i=1:length(xtpt2)
        plot(i,xtpt2{1,i}(1),'r.');
    end
    
    hold off;    
end

    [xtt2 xtpt2 Ptt2 Ptpt2 loglikval2] = kalmanf(kalfilt.sig,kalfilt.F,0,0,kalfilt.Qt,kalfilt.Rt,kalfilt.H);
    plot(kalfilt.sig(1,:),'*');hold on;
    J1=0;
    J2=0;
    for i=1:length(xtpt2)
        plot(i,xtpt2{1,i}(1),'r.');
        J2= J2+trace(inv(Ptt2{1,i}+kalfilt.Qt)*kalfilt.Qt);
        J1 = J1+trace((inv(Ptt2{1,i}+kalfilt.Qt+kalfilt.Rt))*kalfilt.Rt);
    end
    J1 = J1/length(Ptt2);
    J2 = J2/length(Ptt2);
kalfilt.Rt;
kalfilt.Qt;