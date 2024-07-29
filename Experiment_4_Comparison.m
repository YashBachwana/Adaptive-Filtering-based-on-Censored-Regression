clc;

clear all;
close all;


iteration = 1000 ;
g = 0.1 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;
    
    sys_w=[0.1 0.5 0.1]';
    sigma = 1 ;
    sys_tap=zeros(1,3)';
    model_w=zeros(1,3)';
    mu = 0.05 ;
    
    for i=1:length(input)
        sys_tap=[input(i) sys_tap(1:end-1)']';
        if i < 2500
            sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        else 
            sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        end 
        % sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);        
        mdl_opt = cdf('Normal',sys_tap'*model_w,0,sigma) * sys_tap' * model_w + sigma * pdf('Normal',sys_tap'*model_w,0,sigma);
        
        err(i) = sys_opt_cap - mdl_opt;
        
        n_fact = 1 ;
        model_w = model_w + (n_fact * mu * cdf('Normal',sys_tap'*model_w,0,sigma) * sys_tap' * err(i))';
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'b'); hold on ;



mu = 0.01 ; % Update rule coefficient 
N = 3 ;
input_length = 5000 ;
w_sys=[0.1 0.5 0.1] ;
for iter = 1 : 1000

    input = rand(1,5000) - 0.5; % Random signal
    system_noise = awgn(input,30)-input ; % White Gaussian Noise 
    input = [zeros(1,N - 1) input] ; 
    w_LMS = zeros(1,N) ;
    sys_tap = zeros(1,N) ;

    iter
    for i = 1 : input_length 

        sys_tap=[input(i) sys_tap(1:end-1)] ;
        if i < 2500
            sys_opt_cap = cdf('Normal',sys_tap*w_sys',0,1) * sys_tap * w_sys' + pdf('Normal',sys_tap*w_sys',0,1) + system_noise(i);
        else 
            sys_opt_cap = sys_tap * w_sys' + system_noise(i) ;
        end 
        % d = cdf('Normal',sys_tap*w_sys',0,1) * sys_tap * w_sys' + pdf('Normal',sys_tap*w_sys',0,1) + system_noise(i);
        sys_out = w_LMS * sys_tap' ; % System Output
        error = sys_opt_cap - sys_out ; % Error 

        err(i) = error ;

        w_LMS = w_LMS + 2 * mu * error * sys_tap  ; % Update Rule 

    end 
    err_ensemble(iter,:) = err .^ 2 ;
end 

plot(10 * log10(mean(err_ensemble)),'g') ; hold on ; 



















iteration = 1000 ;
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

    mu = 0.5 ;
    
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
        
        model_w1 = model_w1 + (mu * cdf('Normal',sys_tap'*model_w1,0,sigma) * sys_tap' * err1)';
        model_w2 = model_w2 + (mu * err2 * sys_tap) ;
        alpha = alpha + 100 * err(i) * lambda * (1 - lambda) * (y1 - y2) ;
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'r'); hold on ;
legend('CR-LMS','LMS','Convex Combination')

function y = sgmd(x)
    y = 1 / (1 + (exp(-x))) ;
end 