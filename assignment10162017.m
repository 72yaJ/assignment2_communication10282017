% to generate the four symbols by QAM, and compare with the BER
close all;
clear all;
clc;

num = 1e4; % the number of symbols
M = 4; % for 4 PAM
N_bpsym = 2; % for 4 PAM, 2bits/symbol
flag_pic = 09; % set 0 2 4 6 8 10 to show the figure when EbN0 = 0 2 4 6 8 10 
s_s = randi([1 4],1,num)*2-5; % generate the signal of symbols, Am = 2*m-1-M, shoud be -3 -1 1 3
% figure,stem(s_s(1:100))
% return 
t_s = linspace(1,num,num); % time of symbols
% P_d = length(find(s_s==-3))./num % verify P(-3)=P(-1)=P(1)=P(3)=0.25
% P_d = length(find(s_s==-1))./num % verify P(-3)=P(-1)=P(1)=P(3)=0.25
% P_d = length(find(s_s==1))./num  % verify P(-3)=P(-1)=P(1)=P(3)=0.25
% P_d = length(find(s_s==3))./num  % verify P(-3)=P(-1)=P(1)=P(3)=0.25
% return

s_b = zeros(1,N_bpsym*num); % generate the signal of bits
t_b = 1:1/N_bpsym:num+1; % time of bits
t_b = t_b(1:num*N_bpsym);
s_b(find(s_s==-3)*2-1) = -1;
s_b(find(s_s==-3)*2) = -1;
s_b(find(s_s==-1)*2-1) = -1;
s_b(find(s_s==-1)*2) = 1;
s_b(find(s_s==1)*2-1) = 1;
s_b(find(s_s==1)*2) = 1;
s_b(find(s_s==3)*2-1) = 1;
s_b(find(s_s==3)*2) = -1;
t_p(1,:) = ones(1,num)*10; % the pack line
t_p(2,:) = ones(1,num)*-10; % the pack line

BER = zeros(1,7);
SER = zeros(1,7);
s_s1 = zeros(1,length(s_s));
s_b1 = zeros(1,length(s_b));
for m = 1:1:6
    EbN0_db = 2*(m-1);
    Eb = 2.5; % Eb is the average energy per bit, Es=5Ep=5=log2(M)*Eb
    N0_sd = Eb/10.^(EbN0_db/10); % SNR is Eb/N0, N0_sd is noise power spectral density
    N0 = sqrt(N0_sd/2)*randn(1,length(s_s)); % N0_sd/2 is two-sided power spectral density of the noise
    s_sn = s_s+N0; % symbol signal with noise
    s_s1(s_sn<-2) = -3;
    s_s1(s_sn>=-2 & s_sn<0) = -1;
    s_s1(s_sn>=0 & s_sn<2) = 1;
    s_s1(s_sn>=2) = 3;
    SER(m) = length(find(s_s1-s_s))./length(s_s1); % SER of symbol Signal with AWGN

    s_b1(find(s_s1==-3)*2-1) = -1;
    s_b1(find(s_s1==-3)*2) = -1;
    s_b1(find(s_s1==-1)*2-1) = -1;
    s_b1(find(s_s1==-1)*2) = 1;
    s_b1(find(s_s1==1)*2-1) = 1;
    s_b1(find(s_s1==1)*2) = 1;
    s_b1(find(s_s1==3)*2-1) = 1;
    s_b1(find(s_s1==3)*2) = -1;
    BER(m) = length(find(s_b1-s_b))./length(s_b1); % BER of bit Signal with AWGN
    
    if EbN0_db == flag_pic
        figure;
        stem(t_s,s_s,'*'); % Original Symbol-Signal
        hold on;
        stairs(t_b,s_b); % Bit Signal
        hold on;
        stem(t_s,s_sn,'p'); % Symbol Signal+noise
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Original Symbol Signal',...
            'Bit Signal',...
            ['Symbol Signal with ',num2str(EbN0_db),'db AWGN']);
        hold off;
        xlim([1 40]);
        ylim([-4 4]);
        figure;
        stem(t_s,s_s1,'*'); % Symbol Signal with noise after quantization
        hold on;
        stairs(t_b,s_b1); % Bit Signal with noise after quantization
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Symbol Signal with noise after quantization',...
            'Bit Signal with noise after quantization');
        hold off;
        xlim([1 40]);
        ylim([-4 4]);
        figure;
        stem(t_s,s_s,'*'); % Original Symbol Signal
        hold on;
        stem(t_s,s_s1,'o'); % Symbol Signal with noise after quantization
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Original Symbol Signal',...
            'Symbol Signal with noise after quantization');
        hold off;
        xlim([1 40]);
        ylim([-4 4]);
        figure;
        stairs(t_b,s_b,'b'); % Original Bit Signal
        hold on;
        stairs(t_b,s_b1,'r'); % Bit Signal with noise after quantization
        hold on;
        stem(t_s,t_p(1,:),'--y'); % pack line
        hold on;
        stem(t_s,t_p(2,:),'--y'); % pack line
        legend('Original Bit Signal',...
            'Bit Signal with noise after quantization');
        hold off;
        xlim([1 40]);
        ylim([-1.25 1.25]);
    end
end

EbN0_db = linspace(0,10,num);
EbN0 = 10.^(EbN0_db/10);
BER_t_bpsk = qfunc(sqrt(2*EbN0)); % rho = -1,BPSK
% BER_t_ave = 0.5*(1-sqrt(EbN0./(1+EbN0))); % coherence detected BPSK antipodal
% SER_t = (M-1)/M*erfc(sqrt(3/(M^2-1)*log2(M)*EbN0));
[BER_t,SER_t] = berawgn(EbN0_db,'pam',4);
figure;
semilogy(EbN0_db, SER_t,'b');
hold on;
semilogy(EbN0_db, BER_t,'r--');
hold on;
semilogy(EbN0_db, BER_t_bpsk,'r:');
for m = 1:1:m
    hold on;
    plot(2*(m-1),SER(m),'b*');
    hold on;
    plot(2*(m-1),BER(m),'ro');
end
grid on;
ylabel('P');
xlabel('E_b/N_0 (dB)');
legend('theoretical SER in 4PAM','theoretical BER in 4PAM','theoretical BER in BPSK',...
    'simulated SER in 4PAM','simulated BER in 4PAM');
hold off;
return










