%% CIE 428 Project 3 - BFSK
clear;
L=1e3;
Tb=1;
Eb=1;
nc=10;
f1=(nc+1)/(Tb);
f2=(nc+2)/(Tb);
seq=binaryseq(L,1,0);
%% Constellation of the transmitted symbols
hold on
scatter(sqrt(Eb)*seq(seq==1),zeros(1,length(seq(seq==1))),'filled','b');
scatter(zeros(1,length(seq(seq==0))),sqrt(Eb)*(~seq(seq==0)),'filled','b');
title(strcat('Constellation of Transmitted Symbols at E_b =',{' '},'1'));
xlabel('\phi_1');
ylabel('\phi_2');
xlim([0 1.2*sqrt(Eb)]);
ylim([0 1.2*sqrt(Eb)]);
%% Noise Addition
SNR_dB=0;
No=Eb*10^(-SNR_dB/10);
r(1,:)=sqrt(Eb)*seq;
r(2,:)=sqrt(Eb)*(~seq);
r=r+sqrt(No/2)*randn(2,L);
%% Constellation of received symbols
scatter(r(1,:),r(2,:),'filled','r');
title(strcat('Constellation of Received symbols for',{' '},num2str(L),...
{' '},'bits',{' '},'at E_b =',{' '},'1 and SNR =',{' '},num2str(SNR_dB),' dB'));
xlabel('\phi_1');
ylabel('\phi_2');
%% Constellation of transmitted and received symbols
hold on
scatter(sqrt(Eb)*seq(seq==1),zeros(1,length(seq(seq==1))),'filled','b');
scatter(r(1,:),r(2,:),'filled','r');
scatter(zeros(1,length(seq(seq==0))),sqrt(Eb)*(~seq(seq==0)),'filled','b');
title(strcat('Constellation of Transmitted and Received symbols for',{' '},num2str(L),...
{' '},'bits',{' '},'at E_b =',{' '},'1 and SNR =',{' '},num2str(SNR_dB),' dB'));
legend('Transmitted Symbols','Received symbols');
xlabel('\phi_1');
ylabel('\phi_2');
%% BFSK Baseband BER Simulation
L=1e5;
SNR=-3:12;
Eb=1;
ber_sim=zeros(1,length(SNR));
ber=zeros(1,length(SNR));
r=zeros(2,L);
for i=1:length(SNR)
    seq=binaryseq(L,1,0);
    detected=zeros(1,L);
    No=Eb*10^(-SNR(i)/10);
    r(1,:)=sqrt(Eb)*seq;
    r(2,:)=sqrt(Eb)*(~seq);
    r=r+sqrt(No/2)*randn(2,L);
    for j=1:length(seq)
        if(r(1,j)-r(2,j)>0)
            detected(j)=1;
        else
            detected(j)=0;
        end
    end
    ber_sim(i)=1-mean(detected==seq);
    ber(i)=1/2*erfc(sqrt(Eb/(2*No)));
end
hold on;
grid on;
plot(SNR,ber_sim,'-oblue');
plot(SNR,ber,'--diamondred');
ylabel('Bit Error rate of BFSK');
xlabel('Signal-to-Noise ratio in dB');
legend('Simulation','theoretical');
set(gca, 'YScale', 'log');
%% BFSK Modulation - Ensemble Generation
Tb=1;
fs=2*f1*f2;
t=0:1/(fs):Tb;
phi1=sqrt(2/Tb)*cos(2*pi*f1*t);
phi2=sqrt(2/Tb)*cos(2*pi*f2*t);
T=length(t);
N=1e4;
L=30;
X=zeros(N,L*T);
for j=1:N
    seq=binaryseq(L,1,0);
    s_t=zeros(1,T*L);
    for i=1:L
        temp=repelem(seq(i),T);
        s_t((i-1)*T+1:i*T)=((temp).*phi1+(~temp).*phi2);
    end
    X(i,:)=s_t;
end
t=0:1/fs:(length(s_t)-1)/(fs);
plot(t,s_t)
xlabel('t')
ylabel('s_ (t)')
title(strcat('Transmitted FSK Signal for',{' '},num2str(L),' bits at f_1='...
,num2str(f1),', f_2=',num2str(f2),' and Tb=',num2str(Tb)))
%% Baseband PSD of BFSK Signal
Tb=1;
Eb=1;
nc=10;
f1=(nc+1)/(Tb);
f2=(nc+2)/(Tb);
fs=3*f2;
t=0:1/(fs):Tb;
s1=sqrt(2*Eb/Tb)*cos(pi*t/Tb);
s2=sqrt(2*Eb/Tb)*sin(pi*t/Tb);
L=50;
T=length(t);
s_t=zeros(1,T*L);
s1_t=zeros(1,T*L);
N=1e3;
X=zeros(N,T*L);
Y=zeros(N,T*L);
for j=1:N
    seq=binaryseq(L,1,0);
    for i=1:L
        if(seq(i)==1)
            s_t((i-1)*T+1:i*T)=-s2;
            s1_t((i-1)*T+1:i*T)=s1;
        else
            s_t((i-1)*T+1:i*T)=s2;
            s1_t((i-1)*T+1:i*T)=s1;
        end
    end
    X(j,:)=s_t;
    Y(j,:)=s1_t;
end
acf=statacf(X);
acf1=statacf(Y);
psd=abs(fftshift(fft(acf)));
psd1=abs(fftshift(fft(acf1)));
fftsize=length(acf);
f=-fs/2:fs/fftsize:fs/2-fs/fftsize;
psd_total=psd/(2*max(psd))+psd1/max(psd1);
plot(f*Tb,psd_total/max(psd_total),'b');
title('Normalized Power Spectral Density of the transmitted BFSK Signal')
xlabel('fT_b')
ylabel('S(f)')