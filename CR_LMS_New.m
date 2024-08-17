clc;

clear all;
close all;
% System Identification
for itr=1:100

    itr
    input=rand(1,5000)-0.5;
    noise=awgn(input,40)-input;
    
    sys_w=[0.1 0.2 0.1]';

    sigma = 1 ;
    sigma_in = 1 ;
    gamma = zeros(1,3)';
    sys_tap=zeros(1,3)';
    model_w=zeros(1,3)';
    mu = 0.05 ;
    
    for i=1:length(input)
        sys_tap=[input(i) sys_tap(1:end-1)']';

        % sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        if i < 2500
            sys_opt_cap = cdf('Normal',sys_tap'*sys_w,0,1) * sys_tap' * sys_w + pdf('Normal',sys_tap'*sys_w,0,1) + noise(i);
        else 
            sys_opt_cap = sys_tap' * sys_w + noise(i) ;
        end 
        
        mdl_opt = cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * model_w + sigma * pdf('Normal',sys_tap'*gamma,0,sigma);
        
        err(i) = sys_opt_cap - mdl_opt;
        
        if sys_opt_cap > 0 
            lambda = 1 ;
        else 
            lambda = 0 ;
        end 

        if sigma_in >= 1
            slope = 1 ;
        else 
            slope = 0 ; 
        end 
        gamma = gamma + mu * lambda * psi(sys_tap'*gamma,0,sigma) * sys_tap - mu * (1 - lambda) * psi(- sys_tap'*gamma,0,sigma) * sys_tap ;
        
        sigma_in = slope * sigma_in + mu * pdf('Normal',sys_tap'*gamma,0,sigma) * err(i) ;
        sigma = max(sigma_in, 1) ;
        model_w = model_w + (mu * cdf('Normal',sys_tap'*gamma,0,sigma) * sys_tap' * err(i))';
        
    end

err_plot(itr,:)=err.^2;
end
plot(10*log10(mean(err_plot)),'b'); hold on ;

function si = psi(x, mu, sigma)
    % Compute PDF and CDF
    pdf_value = pdf('Normal',x, mu, sigma);
    cdf_value = cdf('Normal',x, mu, sigma);
    
    % Calculate 
    si = pdf_value / cdf_value;
end
