clc ; clear all ; close all ; 

% System Definition 
w_sys = [1,2,1] ;

N = 3 ; % Length of the weight vector 
L = 3 ; % Number of previous vectors

mu = 0.1 ; % Learning Rate
gamma = 0.001 ;
iterations = 100 ; % Iterations 

inp_len = 2000 ;
SNR = 30 ;

for iter = 1 : iterations 

    input = rand(1,inp_len) - 0.5; % Random signal
    % input = [0 input] ;
    system_noise = awgn(input,SNR) - input ; % White Gaussian Noise 

    sys_inp = zeros(L + 1 , N) ;
    curr_inp = zeros(1,N) ;
    w_APA = zeros(1,N) ;

    iter
    for i = 1 : inp_len

        % Creating previous vectors
        curr_inp = [input(i) curr_inp(1:end - 1)] ;
        sys_inp = [curr_inp; sys_inp(1:end - 1, :)];
        
        filter_out = w_APA * sys_inp' ; % Filter Output
        sys_out = w_sys * sys_inp' + system_noise(i) ; % System Output

        error = sys_out - filter_out ; % Error 

        err(i) = w_sys * curr_inp' + system_noise(i) - w_APA * curr_inp' ;

        % w_APA = w_APA + mu * error * sys_inp  ; % Update Rule 
        w_APA = w_APA + mu * error * sys_inp ; % Update Rule 

    end 
    err_ensemble(iter,:) = err .^ 2 ;
end 

plot(10 * log10(mean(err_ensemble)),'r') ; hold on ; 