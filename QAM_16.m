%% CIE 428 Project 3 - 16 QAM
clear;
L=1e3;
E0=1;
seq=sqrt(E0)*quadseq(L,-3,-1,1,3);
%% Constellation of the transmitted symbols
scatter(seq(1,:),seq(2,:),'filled','b');
labels=grayenco(seq,-3*sqrt(E0),-sqrt(E0),sqrt(E0),3*sqrt(E0));
text(seq(1,:),seq(2,:),cellstr(labels),'VerticalAlignment','top','HorizontalAlignment','center');
title(strcat('Constellation of 16-QAM Transmitted Symbols (Gray Encoded) for'...
,{' '},num2str(L),{' '},'Symbols'))
xlabel('\phi_1')
ylabel('\phi_2')
xlim([-4*sqrt(E0) 4*sqrt(E0)])
ylim([-4*sqrt(E0) 4*sqrt(E0)])
%% Noise Addition
SNR_dB=0;
No=E0*10^(-SNR_dB/10);
r=seq+sqrt(No/2)*randn(2,L);
%% Constellation of the Received Symbols
scatter(r(1,:),r(2,:),'filled','r');
title(strcat('Constellation of',{' '},num2str(L),{' '},'Received 16-QAM Symbols',{' '},...
'at E_0 =',{' '},num2str(E0),{' '},'and SNR =',{' '},num2str(SNR_dB),{' '},'dB'))
xlabel('\phi_1')
ylabel('\phi_2')
%% Constellation of transmitted and received symbols
hold on
scatter(seq(1,:),seq(2,:),'filled','b');
scatter(r(1,:),r(2,:),'filled','r');
xlabel('\phi_1')
ylabel('\phi_2')
legend('Transmitted Symbols','Received Symbols')
title(strcat('Constellation of',{' '},num2str(L),{' '},'Transmitted and Received 16-QAM Symbols',{' '},...
'at E_0 =',{' '},num2str(E0),{' '},'and SNR =',{' '},num2str(SNR_dB),{' '},'dB'))
%% Baseband simulation for BER
clear;
L=1e5;
SNR=-3:12;
E0=1;
M=16;
Eav=2/3*(M-1)*E0;
ber_sim=zeros(1,length(SNR));
ber=zeros(1,length(SNR));
symbols=sqrt(E0)*[-3 -3 -3 -3 -1 -1 -1 -1 1 1 1 1 3 3 3 3;...
-3 -1 1 3 -3 -1 1 3 -3 -1 1 3 -3 -1 1 3];
for i=1:length(SNR)
    seq=sqrt(E0)*quadseq(L,-3,-1,1,3);
    detected=zeros(2,L);
    No=1/4*Eav*10^(-SNR(i)/10);
    r=seq+sqrt(No/2)*randn(2,L);
    for j=1:length(seq)
        distances=zeros(1,length(symbols));
        for k=1:length(symbols)
            distances(k)=norm(r(:,j)-symbols(:,k));
        end
        [~,I]=min(distances);
        detected(:,j)=symbols(:,I);
    end
    temp=int8(detected(1,:)==seq(1,:))+int8(detected(2,:)==seq(2,:));
    ber_sim(i)=(1-mean(temp==2));
    ber(i)=3/2*erfc(sqrt(1/10*Eav/No))-9/16*(erfc(sqrt(1/10*Eav/No)))^2;
end
hold on;
grid on;
plot(SNR,ber_sim,'-oblue');
plot(SNR,ber,'--diamondred');
ylabel('Bit Error rate of QAM');
xlabel('Signal-to-Noise ratio in dB');
legend('Simulation','theoretical');
set(gca, 'YScale', 'log');
%% 16 - QAM Modulation - Ensemble Generation
L=10;
N=1e3;
fc=1e2;
fs=4*fc;
T=0.5;
E0=1;
t=0:1/fs:T;
phi1=sqrt(2/T)*cos(2*pi*fc*t);
phi2=sqrt(2/T)*sin(2*pi*fc*t);
X=zeros(N,length(t)*L);
for i=1:N
    seq=quadseq(L,-3,-1,1,3);
    s_t=zeros(1,length(t)*L);
    for j=1:L
        s_t((j-1)*length(t)+1:j*length(t))=s_t((j-1)*length(t)+1:j*length(t))+...
    sqrt(E0)*seq(1,j)*phi1-sqrt(E0)*seq(2,j)*phi2;
    end
    X(i,:)=s_t;
end
plot(0:1/fs:(length(s_t)/fs)-1/fs,s_t);
title(strcat('16-QAM Modulated Signal at fc=',{' '},num2str(fc),{' '}...
,'and E_o=',{' '},num2str(E0),' for',{' '},num2str(L),{' '},'bits'))
xlabel('t')
ylabel('s(t)')
%% PSD of Passband Signal
acf=statacf(X);
psd=abs(fftshift(fft(acf)));
f=-fs/2:fs/length(psd):fs/2-1/length(psd);
plot(f,psd/max(psd));
title('Normalized PSD of Passband 16-QAM Signal on a logarithmic scale')
ylabel('Normalized Power Spectral Density')
xlabel('Frequency (Hz)')
set(gca, 'YScale', 'log');
%% Baseband PSD of 16-QAM Signal
L=1e3;
N=1e4;
T=5;
E0=1;
X=zeros(N,L*T);
for i=1:2:N-1
    seq=sqrt(E0)*quadseq(L,-3,-1,1,3);
    X(i,:)=repelem(seq(1,:),T);
    X(i+1,:)=repelem(seq(2,:),T);
end
acf=statacf(X);
psd=abs(fftshift(fft(acf)));
f=-fs/2:fs/length(psd):fs/2-1/length(psd);
plot(f,psd/max(psd));
ylabel('Normalized PSD of 16-QAM Signal')
xlabel('Frequency (Hz)')