clc ; clear all ; close all ; 

iteration = 1000 ;
g = 0.1 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;
    
    sys_w=[0.1 0.5 0.1]';
    sigma = 1 ;
    lambda = 0.5 ;
    alpha = 0 ;
    sys_tap=zeros(1,3)';

    model_w1=zeros(1,3)';
    model_w2=zeros(1,3)';

    mu = 0.1 ;
    
    for i=1:length(input)
        sys_tap=[input(i) sys_tap(1:end-1)']';
        if i < 2500
            sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        else 
            sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        end 
        % sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        % sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        
        lambda = sgmd(alpha) ;

        y1 = cdf('Normal',sys_tap'*model_w1,0,sigma) * sys_tap' * model_w1 + sigma * pdf('Normal',sys_tap'*model_w1,0,sigma);
        y2 = sys_tap'*model_w2 ;
        
        mdl_opt = (lambda * y1) + ((1 - lambda) * y2) ;

        err1 = sys_opt_cap - y1 ;
        err2 = sys_opt_cap - y2 ;
        err(i) = sys_opt_cap - mdl_opt;
        
        n_fact1 = 1  / (((cdf('Normal',sys_tap'*model_w1,0,sigma) ^ 2) * (sys_tap' * sys_tap)) + g) ;
        n_fact2 = 1 / ((sys_tap' * sys_tap) + g) ;
        model_w1 = model_w1 + n_fact1 * (mu * cdf('Normal',sys_tap'*model_w1,0,sigma) * sys_tap' * err1)';
        model_w2 = model_w2 + n_fact2 * (mu * err2 * sys_tap) ;
        alpha = alpha + err(i) * lambda * (1 - lambda) * (y1 - y2) ;
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'r'); hold on ;












clear ;

iteration = 1000 ;
g = 0.1 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;
    
    sys_w=[0.1 0.5 0.1]';
    sigma = 1 ;
    lambda = 0.5 ;
    alpha = 0 ;
    sys_tap=zeros(1,3)';

    model_w1=zeros(1,3)';
    model_w2=zeros(1,3)';

    mu = 0.1 ;
    
    for i=1:length(input)
        sys_tap=[input(i) sys_tap(1:end-1)']';
        if i < 2500
            sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        else 
            sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        end 
        % sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        % sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        
        lambda = sgmd(alpha) ;

        y1 = cdf('Normal',sys_tap'*model_w1,0,sigma) * sys_tap' * model_w1 + sigma * pdf('Normal',sys_tap'*model_w1,0,sigma);
        y2 = sys_tap'*model_w2 ;
        
        mdl_opt = (lambda * y1) + ((1 - lambda) * y2) ;

        err1 = sys_opt_cap - y1 ;
        err2 = sys_opt_cap - y2 ;
        err(i) = sys_opt_cap - mdl_opt;
        
        n_fact1 = 1 ;
        n_fact2 = 1 ; 
        model_w1 = model_w1 + n_fact1 * (mu * cdf('Normal',sys_tap'*model_w1,0,sigma) * sys_tap' * err1)';
        model_w2 = model_w2 + n_fact2 * (mu * err2 * sys_tap) ;
        alpha = alpha + err(i) * lambda * (1 - lambda) * (y1 - y2) ;
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'b'); hold on ;




function y = sgmd(x)
    y = 1 / (1 + (exp(-x))) ;
end 