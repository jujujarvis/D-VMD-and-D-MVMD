function [u, u_hat, omega] = D_MVMD(signal, alpha, tau1, tau2, K, DC, init, tol)
% De-mixing (multivariate) variational mode decomposition (D-VMD, D-MVMD)
% Author: Jiawei Jian 
% This function is modified from Variational Mode Decomposition 
% by Konstantin Dragomiretskiy and Dominique Zosso.

%% Citation
% When using this code, please do cite our paper: 
% J. Jian, Z. Lu, J. Liu, et al. "De-Mixing Variational Mode Decomposition 
% and its Application on Operational Modal Analysis 
% in the Presence of Closely Spaced Modes", Measurement. 

%% Code description
% Input and Parameters:
% ---------------------
% signal  - the time domain multivariate signal to be decomposed
% alpha   - the penalty factor as well as the bandwidth control parameter
% tau1    - the ascent parameter for first Lagrangian multiplier 
%           (i.e., noise tolerance parameter)
% tau2    - the ascent parameter for second Lagrangian multiplier 
%           (i.e., correlation tolerance parameter)
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%           1 = all omegas start uniformly distributed
%           2 = all omegas initialized randomly
% tol     - tolerance of convergence criterion; typically around 1e-6
%
% Output:
% -------
% u       - the decomposed modes
% u_hat   - the spectra of the modes
% omega   - estimated mode center-frequencies

%% Preparations
% signal length, channel number, sampling frequency
save_T = size(signal,1); 
CC = size(signal,2); 
fs = 1/save_T; 

% extend the signal by mirroring
T = save_T;
f_mirror(1:T/2,:) = signal(T/2:-1:1,:);
f_mirror(T/2+1:3*T/2,:) = signal;
f_mirror(3*T/2+1:2*T,:) = signal(T:-1:T/2+1,:);
f = f_mirror;

% time domain 0 to T of mirrored signal
T = size(f,1);
t = (1:T)/T;

% spectral domain 
freqs = t-0.5-1/T;

% maximum number of iterations
N = 500;

% for future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat.';
f_hat_plus(:,1:T/2) = 0;

% matrix keeping track of every iterant of spectra
u_hat_plus = rand(CC, T, K);
u_hat_plus(:, 1:T/2, :) = 0;

% initialization of omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end

% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(CC, length(freqs));
lambda2 = zeros(CC, K);

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = zeros(CC,T,K); % accumulator

%% Main loop for iterative updates

while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    u_hat_plus_last=u_hat_plus;
    lambda_hat_last=lambda_hat;
    lambda2_last=lambda2;

    % update first mode accumulator
    k = 1;
    for c=1:CC
        % accumulator
        sum_uk(c,:,k) = sum(u_hat_plus(c,:,:),3) - u_hat_plus_last(c,:,k);
        
        % update three scalar items
        sca1=lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
            sum( f_hat_plus(c,T/2+1:T)-2*u_hat_plus(c,T/2+1:T,k) );
        sca2=-lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))^3/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
            ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' )*sum(u_hat_plus(c,T/2+1:T,k));
        sca3=lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))^3*...
            ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' )*sum(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k));
        
        % update spectrum of the first mode
        u_hat_plus(c,:,k) = ( f_hat_plus(c,:) - sum_uk(c,:,k) - lambda_hat(c,:)/2 - ...
            sca1-sca2-sca3 )./ (1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);   

        % update lambda2
        lambda2(c,k)=lambda2_last(c,k)+tau2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
            ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' );
    end

    % update first omega
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*sum((abs(u_hat_plus(:, T/2+1:T, k)).^2),1).')/sum(sum(abs(u_hat_plus(:,T/2+1:T,k)).^2));
    end
    
    % update of any other mode
    for k=2:K
        for c=1:CC
            % accumulator
            sum_uk(c,:,k) = sum(u_hat_plus(c,:,:),3) - u_hat_plus_last(c,:,k);
            
            % update three scalar items
            sca1=lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
                sum( f_hat_plus(c,T/2+1:T)-2*u_hat_plus(c,T/2+1:T,k) );
            sca2=-lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))^3/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
                ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' )*sum(u_hat_plus(c,T/2+1:T,k));
            sca3=lambda2(c,k)/2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))^3*...
                ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' )*sum(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k));
            
            % update spectra
            u_hat_plus(c,:,k) = ( f_hat_plus(c,:) - sum_uk(c,:,k) - lambda_hat(c,:)/2 - ...
                sca1-sca2-sca3 )./ (1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);   
            
            % update lambda2
            lambda2(c,k)=lambda2_last(c,k)+tau2/norm(u_hat_plus(c,T/2+1:T,k))/norm(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))*...
                ( u_hat_plus(c,T/2+1:T,k)*(f_hat_plus(c,T/2+1:T)-u_hat_plus(c,T/2+1:T,k))' );
        end        

        % update center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*sum((abs(u_hat_plus(:, T/2+1:T, k)).^2),1).')/sum(sum(abs(u_hat_plus(:,T/2+1:T,k)).^2));
        
    end
    
    % update lambda1
    for c=1:CC
        lambda_hat(c,:) = lambda_hat_last(c,:) + tau1*(sum(u_hat_plus(c,:,:),3) - f_hat_plus(c,:));
    end

    % loop counter
    n = n+1;
    
    % convergence check
    uDiff = eps;
    for i=1:K
        for c=1:CC
            uDiff = uDiff + 1/T*(u_hat_plus(c,:,i)-u_hat_plus_last(c,:,i))*conj((u_hat_plus(c,:,i)-u_hat_plus_last(c,:,i)))';
        end
    end
    uDiff = abs(uDiff);
    
end

%% Postprocessing and cleanup

% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(CC,T, K);
u_hat(:,(T/2+1):T,:) = u_hat_plus(:,(T/2+1):T,:);
u_hat(:,(T/2+1):-1:2,:) = conj(u_hat_plus(:,(T/2+1):T,:));
u_hat(:,1,:) = conj(u_hat(:,end,:));

u = zeros(K,length(t),CC);
for k = 1:K
    for c=1:CC
        u(k,:,c)=real(ifft(ifftshift(u_hat(c,:,k))));
    end
end

% remove mirror part
u = u(:,T/4+1:3*T/4,[floor(CC/2)+1:CC,1:floor(CC/2)]);

% recompute spectra
clear u_hat;
u_hat = zeros(K,size(u,2),CC);
for k = 1:K
    for c=1:CC
        u_hat(k,:,c)=fftshift(fft(u(k,:,c)));
    end
end

end
