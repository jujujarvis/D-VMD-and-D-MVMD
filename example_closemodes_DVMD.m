clc
clear all

%% simulating signal with closely spaced modes
SampleFrequency=1000;
dt=0.001;
tt=dt:dt:15;
N=length(tt);
df=1/tt(end);
ff=(0:N-1)*df;

th_Freq=[5 15 15.2];
s1=exp(-th_Freq(1)*0.005*2*pi*tt).*sin(th_Freq(1)*2*pi*tt);
s2=exp(-th_Freq(2)*0.005*2*pi*tt).*sin(th_Freq(2)*2*pi*tt);
s3=exp(-th_Freq(3)*0.005*2*pi*tt).*sin(th_Freq(3)*2*pi*tt);

B1=0.005*2*th_Freq(2)*2*pi;
B2=0.005*2*th_Freq(3)*2*pi;
B=(B1+B2)/2;
MOF=B/abs(th_Freq(3)-th_Freq(2));

A=[0.8 -0.6 0.7;0.5 0.7 -0.4;0.3 0.5 0.6];
ss=(s1.*A(:,1)+s2.*A(:,2)+s3.*A(:,3))';
signal=ss;

figure
subplot(1,2,1)
for i=1:size(signal,2)
    plot(tt,signal(:,i));hold on;xlim([0 15]);%ylim([-10 10]);
end

subplot(1,2,2)
for i=1:size(signal,2)
    plot(ff,abs(fft(signal(:,i))));hold on;xlim([0 20]);%ylim([-10 10]);
end

%% parameters setting
alpha = 50000;        % the bandwidth control parameter
tau1 = 0.01;          % the ascent parameter 1
tau2 = 0.01;          % the ascent parameter 2
KK = 3;               % mode number
DC = 0;               % DC mode parameter
init = 1;             % initialize omegas
tol = 1e-6;           % convergence criterion

%% mode decomposition
tic
[u, u_hat, omega] = D_MVMD(signal, alpha, tau1, tau2, KK, DC, init, tol);
omega=omega*SampleFrequency;
toc

tic
[u2, u_hat2, omega2] = MVMD(signal, alpha, tau1, KK, DC, init, tol);
omega2=omega2*SampleFrequency;
toc

num_omega=size(omega,1);
num_omega2=size(omega2,1);
[~,order1] = sort(omega(num_omega,:));
omega=omega(:,order1);
[~,order2] = sort(omega2(num_omega2,:));
omega2=omega2(:,order2);

%% drawing
col=['b','r','g'];

figure
for num_m=1:size(u,1)
    subplot(size(u,1),2,2*num_m-1)
    for i=1:size(u,3)
        sep_signal=u(order1(num_m),:,i);
        sep_signal=sep_signal';
        plot(tt,sep_signal,col(i));hold on;xlim([0,10])
    end
    subplot(size(u,1),2,2*num_m)
    for i=1:size(u,3)
        sep_signal=u(order1(num_m),:,i);
        sep_signal=sep_signal';
        plot(ff,abs(fft(sep_signal)),col(i));hold on;xlim([0,20])
    end
    legend('Channel 1','Channel 2','Channel 3');
end

figure
for num_m=1:size(u2,1)
    subplot(size(u2,1),2,2*num_m-1)
    for i=1:size(u2,3)
        sep_signal=u2(order2(num_m),:,i);
        sep_signal=sep_signal';
        plot(tt,sep_signal,col(i));hold on;xlim([0,10])
    end
    subplot(size(u2,1),2,2*num_m)
    for i=1:size(u2,3)
        sep_signal=u2(order2(num_m),:,i);
        sep_signal=sep_signal';
        plot(ff,abs(fft(sep_signal)),col(i));hold on;xlim([0,20])
    end
    legend('Channel 1','Channel 2','Channel 3');
end
