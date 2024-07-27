clc;

clear all;
close all;
iteration = 100 ;
g = 0.1 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;
    
    sys_w=[0.1 0.5 0.1]';
    sigma = 1 ;
    % gamma = zeros(1,3)'; %____________________________________________
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
        
        n_fact = 1  / (((cdf('Normal',sys_tap'*model_w,0,sigma) ^ 2) * (sys_tap' * sys_tap)) + g) ;
        model_w = model_w + (n_fact * mu * cdf('Normal',sys_tap'*model_w,0,sigma) * sys_tap' * err(i))';
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'r'); hold on ;






iteration = 100 ;
g = 0.1 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;
    
    sys_w=[0.1 0.5 0.1]';
    sigma = 1 ;
    % gamma = zeros(1,3)'; %____________________________________________
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





clear ;


mu = 0.01 ; % Update rule coefficient 
N = 3 ;
input_length = 5000 ;
w_sys=[0.1 0.5 0.1] ;
for iter = 1 : 100

    input = rand(1,5000) - 0.5; % Random signal
    system_noise = awgn(input,35)-input ; % White Gaussian Noise 
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


legend('CR-NLMS','CR-LMS','LMS')