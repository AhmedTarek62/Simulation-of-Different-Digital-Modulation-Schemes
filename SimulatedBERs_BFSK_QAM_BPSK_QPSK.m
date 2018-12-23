%% BER of BFSK - 16QAM - BPSK - QPSK
% BFSK Baseband BER Simulation
clear;
L=1e5;
SNR=-3:12;
Eb=1;
ber_bfsk=zeros(1,length(SNR));
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
    ber_bfsk(i)=1-mean(detected==seq);
end
% 16-QAM Baseband simulation for BER
L=1e5;
SNR=-3:12;
E0=1;
M=16;
Eav=2/3*(M-1)*E0;
ber_qam=zeros(1,length(SNR));
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
    ber_qam(i)=(1-mean(temp==2));
end
% BPSK Baseband simulation for BER
SNR=-3:1:12;
L=1e6;
Ep=1;
seq=binaryseq(L,0,1);
encoded=PolarNRZ(seq,1)*sqrt(Ep);
ber_bpsk=zeros(1,length(SNR));
for i=1:length(SNR)
    var=Ep*10^(-SNR(i)/10);
    x_t=encoded+ sqrt(var/2)*randn(1,L);
    detected=x_t>=0;
    ber_bpsk(i)=1-mean(detected==seq);
end
% QPSK Baseband simulation for BER
SNR=-3:1:12;
L=1e6;
Ep=1;
seq=binaryseq(L,0,1);
seq1=binaryseq(L,0,1);
encoded=PolarNRZ(seq,1)*sqrt(Ep/2);
encoded1=PolarNRZ(seq1,1)*sqrt(Ep/2);
ber_qpsk=zeros(1,length(SNR));
for i=1:length(SNR)
    var=(Ep/2)*10^(-SNR(i)/10);
    x_t=encoded+ sqrt(var/2)*randn(1,L);
    x_t1=encoded1+ sqrt(var/2)*randn(1,L);
    detected=x_t>=0;
    detected1=x_t1>=0;
    errors=0;
    for j=1:L
        if detected(j)~=seq(j) && detected1(j)~=seq1(j)
            errors=errors+2;
        elseif detected(j)~=seq(j) || detected1(j)~=seq1(j)
            errors=errors+1;
        end
    end
    ber_qpsk(i)=errors/(2*length(detected));
end
%% 
hold on
grid on
plot(SNR,ber_bpsk,'r-+');
plot(SNR,ber_qpsk,'b-*');
plot(SNR,ber_bfsk,'k-diamond');
plot(SNR,ber_qam,'m-x');
legend('BPSK','QPSK','BFSK','16-QAM')
title('Bit Error Rate simulation of Different Digital Modulation Schemes')
ylabel('Bit Error rate');
xlabel('Signal-to-Noise ratio in dB');
set(gca, 'YScale', 'log');

