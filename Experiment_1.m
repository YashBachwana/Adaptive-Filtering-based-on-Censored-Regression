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
    gamma = zeros(1,3)';
    sys_tap=zeros(1,3)';
    model_w=zeros(1,3)';
    mu = 0.05 ;
    
    for i=1:length(input)
        if i == 2500
            sys_w = [0.5 0.4 0.3]' ;
        else 

        end 
        sys_tap=[input(i) sys_tap(1:end-1)']';

        sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        % if sys_opt > 0
        %     sys_opt_cap = sys_opt ; 
        % else 
        %     sys_opt_cap = 0 ;
        % end 
        
        mdl_opt = cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * model_w + sigma * pdf('Normal',sys_tap'*gamma,0,sigma);
        
        err(i) = sys_opt_cap - mdl_opt;
        
        if sys_opt_cap > 0 
            lambda = 1 ;
        else 
            lambda = 0 ;
        end 
        
        n_fact = 1  / (((cdf('Normal',sys_tap'*model_w,0,sigma) ^ 2) * (sys_tap' * sys_tap)) + g) ;
        gamma = gamma + n_fact * mu * lambda * psi(sys_tap'*gamma,0,sigma) * sys_tap - mu * (1 - lambda) * psi(- sys_tap'*gamma,0,sigma) * sys_tap ;
        sigma = sigma + n_fact * mu * pdf('Normal',sys_tap'*gamma,0,sigma) * err(i) ;
        model_w = model_w + (n_fact * mu * cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * err(i))';
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'r'); hold on ;












clear;
iteration = 100 ;
% System Identification
for itr=1:iteration

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,30)-input;

    sys_w=[0.1 0.5 0.1]';

    sigma = 1 ;
    gamma = zeros(1,3)';
    sys_tap=zeros(1,3)';
    model_w=zeros(1,3)';
    mu = 0.05 ;

    for i=1:length(input)
        if i == 2500
            sys_w = [0.5 0.4 0.3]' ;
        else 

        end 
        sys_tap=[input(i) sys_tap(1:end-1)']';

        sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        % if sys_opt > 0
        %     sys_opt_cap = sys_opt ; 
        % else 
        %     sys_opt_cap = 0 ;
        % end 

        mdl_opt = cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * model_w + sigma * pdf('Normal',sys_tap'*gamma,0,sigma);

        err(i) = sys_opt_cap - mdl_opt;

        if sys_opt_cap > 0 
            lambda = 1 ;
        else 
            lambda = 0 ;
        end 

        gamma = gamma + mu * lambda * psi(sys_tap'*gamma,0,sigma) * sys_tap - mu * (1 - lambda) * psi(- sys_tap'*gamma,0,sigma) * sys_tap ;
        sigma = sigma + mu * pdf('Normal',sys_tap'*gamma,0,sigma) * err(i) ;
        model_w = model_w + (mu * cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * err(i))';
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
    X = zeros(1,N) ;

    iter
    for i = 1 : input_length 
            
        if i == 2500
            sys_w = [0.5 0.4 0.3] ;
        else 

        end 
        X = [input(i) X(1:end - 1)] ;
        d = cdf('Normal',X*w_sys',0,1) * X * w_sys' + pdf('Normal',X*w_sys',0,1) + system_noise(i);
        sys_out = w_LMS * X' ; % System Output
        error = d - sys_out ; % Error 

        err(i) = error ;

        w_LMS = w_LMS + 2 * mu * error * X  ; % Update Rule 

    end 
    err_ensemble(iter,:) = err .^ 2 ;
end 

plot(10 * log10(mean(err_ensemble)),'g') ; hold on ; 
legend('CR-NLMS','CR-LMS','LMS')








function si = psi(x, mu, sigma)
    % Compute PDF and CDF
    pdf_value = pdf('Normal',x, mu, sigma);
    cdf_value = cdf('Normal',x, mu, sigma);
    
    % Calculate 
    si = pdf_value / cdf_value;
end