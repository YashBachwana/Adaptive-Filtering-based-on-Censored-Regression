clc ; clear ; close all ;

mu = 0.1 ; % Learning Rate 
iterations = 50 ; % Iterations 
w_sys = [1,2,1] ;

inp_len = 1000 ;
SNR = 30 ;
N = 3 ; 

% Kernel function (Gaussian kernel)
sigma = 10; % Kernel width
kernel = @(x, y) exp(-norm(x - y)^2 / (2 * sigma^2));



for iter = 1 : iterations
    input = rand(1, inp_len) - 0.5;
    Centers = []; 
    Coefficients = [];
    % Generate system noise using AWGN
    system_noise = awgn(input, SNR) - input;
    sys_inp = zeros(1,N) ;
    iter
    for i = 1:inp_len
        sys_inp = [input(i) sys_inp(1:end - 1)] ;

        if isempty(Centers)
            y = 0;
        else
            k = arrayfun(@(j) kernel(sys_inp, Centers(j)), 1:length(Centers));
            y = sum(Coefficients .* k);
        end
    
        % Calculate the error
        sys_out = w_sys * sys_inp' + system_noise(i) ;

        % error = exp(-norm(sys_inp)^2 / (2 * sigma^2)) + system_noise(i) - y;
        error = sys_out - y ;
        err(i) = error;
    
        Centers = [Centers, input(i)];
        Coefficients = [Coefficients, mu * error];        
    end 

    err_ensemble(iter,:) = err .^ 2 ;

end 

plot(10 * log10(mean(err_ensemble)),'r') ; hold on ; 