clc ; clear all ; close all ; 

% System Definition 
w_sys = [1,20,1,30] ;

N = 4 ; % Length of the weight vector 


mu = 0.1 ; % Learning Rate 
gamma = 0.1 ;
iterations = 100 ; % Iterations 

inp_len = 2000 ;
SNR = 30 ;

for iter = 1 : iterations 

    input = rand(1,inp_len) - 0.5; % Random signal
    system_noise = awgn(input,SNR) - input ; % White Gaussian Noise 

    sys_inp = zeros(1,N) ;
    w_NLMS = zeros(1,N) ;

    iter
    for i = 1 : inp_len
        sys_inp = [input(i) sys_inp(1:end - 1)] ;

        filter_out = w_NLMS * sys_inp' ; % Filter Output
        sys_out = w_sys * sys_inp' + system_noise(i) ; % System Output

        error = sys_out - filter_out ; % Error 

        err(i) = error ;
        n_fact = 1 / (gamma + (sys_inp * sys_inp')) ;
        w_NLMS = w_NLMS + n_fact * (mu) * (error * sys_inp)  ; % Update Rule 

    end 
    err_ensemble(iter,:) = err .^ 2 ;
end 

plot(10 * log10(mean(err_ensemble)),'r') ; hold on ; 
