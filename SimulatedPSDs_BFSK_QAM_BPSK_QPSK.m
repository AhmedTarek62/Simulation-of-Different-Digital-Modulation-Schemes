% Baseband PSD of BFSK Signal
clear;
Tb=2;
Eb=1;
nc=0;
f1=(nc+1)/(Tb);
f2=(nc+2)/(Tb);
fs=2*(f1+f2);
t=0:1/(fs):Tb;
s1=sqrt(2*Eb/Tb)*cos(pi*t/Tb);
s2=sqrt(2*Eb/Tb)*sin(pi*t/Tb);
L=1e3;
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
psd=abs(fftshift(fft(acf)))/(2*Eb);
psd1=abs(fftshift(fft(acf1)))/(2*Eb);
fftsize=length(acf);
f_bfsk=(-fs/2:fs/fftsize:fs/2-fs/fftsize)*(Tb);
% PSD of Baseband 16-QAM Signal
L=1e3;
N=1e3;
T=6;
fs=1/T;
E0=1;
Eav=(2/3)*E0*(15);
X=zeros(N,L*T);
for i=1:2:N-1
    seq=sqrt(E0)*quadseq(L,-3,-1,1,3);
    X(i,:)=repelem(seq(1,:),T);
    X(i+1,:)=repelem(seq(2,:),T);
end
acf_qam=statacf(X);
psd_qam=abs(fftshift(fft(acf_qam)))/(Eav/3);
f_qam=(-fs/2:fs/length(psd_qam):fs/2-fs/length(psd_qam))*T/4;
% PSD of BPSK Baseband Signal
N=1e3;
L=1e3;
Tp=3;
fs=1/Tp;
Ep=1;
s=zeros(N,L*Tp);
for i=1:N
    seq=binaryseq(L,0,1);
    encoded=sqrt(Ep)*PolarNRZ(seq,Tp);%PNRZ Signal
    s(i,:)=encoded;
end
acf_bpsk=statacf(s);
f_bpsk=(-fs/2:fs/length(acf_bpsk):fs/2-fs/length(acf_bpsk))*Tp;
psd_bpsk=abs(fftshift(fft(acf_bpsk)))/(2*Ep);
% PSD of Baseband QPSK Signal
N=1e3;
L=1e3;
T=6;
fs=1/T;
E=2;
s=zeros(N,L*T);
for i=1:N
    seq=binaryseq(L,0,1);
    encoded=sqrt(E/2)*PolarNRZ(seq,T);%PNRZ Signal
    s(i,:)=encoded;
end
acf_qpsk=statacf(s);
f_qpsk=(-fs/2:fs/length(acf_qpsk):fs/2-fs/length(acf_qpsk))*T/2;
psd_qpsk=abs(fftshift(fft(acf_qpsk)))/(E);
%% 
maximum=max(horzcat(psd_qam,psd_bpsk,psd_qpsk,psd));
hold on
plot(f_bfsk,psd/(2.5*maximum)+psd1/max(psd1),'m');
plot(f_qam,psd_qam/maximum,'r');
plot(f_bpsk,psd_bpsk/maximum,'b')
plot(f_qpsk,psd_qpsk/maximum,'y');
legend('BFSK','16-QAM','BPSK','QPSK');
title('Normalized PSD for Different Digital Modulation Schemes')
ylabel('Power Spectral Density')
xlabel('Frequency')
