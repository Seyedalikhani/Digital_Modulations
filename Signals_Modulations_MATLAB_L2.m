clc
clear
close all
% signal generation example
% EITG05, Michael Lentmaier, September 2020

text='digital';
m=double(text); % ASCII integer representation
b=reshape(de2bi(m,7).',1,[]); % binary sequence, 7 bit per character 

b(50)=0;
b(51)=0;
b(52)=0;




%b=[1 1 1 1];
%b=randi(2,1,500)-1;

fsamp=100; % number samples per Ts
Fcarrier=4*fsamp;  % Carrier Frequency
Eg=1/fsamp;        % Eg (Energy per Bit)

t=0:1/fsamp:length(b)-1/fsamp;
t_samp=0:1/fsamp:1-1/fsamp;






% (((((((((((((((((( BPSK ))))))))))))))))))))))

BPSK_Baseband=[];
BPSK_Modulated=[];

for k=1:length(b)
    if (b(k)==0)
        S=-ones(1,fsamp);
        S_M=S.*cos(2*pi*Fcarrier*t_samp/fsamp);
    else
        S=ones(1,fsamp);
        S_M=S.*cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    BPSK_Baseband=[BPSK_Baseband S];
    BPSK_Modulated=[BPSK_Modulated S_M];
end


figure(1)
ax1 =subplot(2,1,1);
plot(t,BPSK_Baseband,'linewidth',1.5);
ylim([-2 2])
xlabel('$t/T_s$','interpreter','latex');
ylabel('$s(t)$','interpreter','latex');
title('BPSK BaseBand Signal')
grid on;
ax2 =subplot(2,1,2);
plot(t,BPSK_Modulated,'linewidth',1.5);
title('BPSK Modulated Signal')
ylim([-2 2])
linkaxes([ax1,ax2],'x');
grid on;





% (((((((((((((((((( QPSK )))))))))))))))))))))

QPSK_Modulated=[];

t_samp=0:1/fsamp:2-1/fsamp;

for k=1:2:length(b)
    if (b(k)==0 && b(k+1)==0)
        S_M=cos(2*pi*Fcarrier*t_samp/fsamp+2*pi*0/4);
    end
    if (b(k)==1 && b(k+1)==0)
        S_M=cos(2*pi*Fcarrier*t_samp/fsamp+2*pi*1/4);
    end
    if (b(k)==0 && b(k+1)==1)
        S_M=cos(2*pi*Fcarrier*t_samp/fsamp+2*pi*2/4);
    end
    if (b(k)==1 && b(k+1)==1)
        S_M=cos(2*pi*Fcarrier*t_samp/fsamp+2*pi*3/4);
    end
    QPSK_Modulated=[QPSK_Modulated S_M];
end

figure(2)
ax1 =subplot(2,1,1);
plot(t,BPSK_Baseband,'linewidth',1.5);
ylim([-2 2])
xlabel('$t/T_s$','interpreter','latex');
ylabel('$s(t)$','interpreter','latex');
title('BPSK BaseBand Signal')
grid on;
ax2 =subplot(2,1,2);
plot(t,QPSK_Modulated,'linewidth',1.5);
title('QPSK Modulated Signal')
ylim([-2 2])
linkaxes([ax1,ax2],'x');
grid on;



% Average Energy of QPSK - Evaluation
M=4;
fun = @(t) (cos(2*pi*Fcarrier*t/fsamp+2*pi*0/4)).^2;
q1 = integral(fun,0,1/fsamp);
fun = @(t) (cos(2*pi*Fcarrier*t/fsamp+2*pi*1/4)).^2;
q2 = integral(fun,0,1/fsamp);
fun = @(t) (cos(2*pi*Fcarrier*t/fsamp+2*pi*2/4)).^2;
q3 = integral(fun,0,1/fsamp);
fun = @(t) (cos(2*pi*Fcarrier*t/fsamp+2*pi*3/4)).^2;
q4 = integral(fun,0,1/fsamp);

Es_avg_QPSK_Evaluation=1/M*(q1+q2+q3+q4)


% Average Energy of QPSK - Formula
Es_avg_QPSK_Formula=Eg/2


disp('****************************************************************')


% (((((((((((((((((( 4-PAM )))))))))))))))))))))

PAM_Modulated=[];

t_samp=0:1/fsamp:2-1/fsamp;

for k=1:2:length(b)
    if (b(k)==0 && b(k+1)==0)
        S_M=-3*cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==1 && b(k+1)==0)
        S_M=-1*cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==0 && b(k+1)==1)
        S_M=1*cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==1 && b(k+1)==1)
        S_M=3*cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    PAM_Modulated=[PAM_Modulated S_M];
end

figure(3)
ax1 =subplot(2,1,1);
plot(t,BPSK_Baseband,'linewidth',1.5);
ylim([-2 2])
xlabel('$t/T_s$','interpreter','latex');
ylabel('$s(t)$','interpreter','latex');
title('BPSK BaseBand Signal')
grid on;
ax2 =subplot(2,1,2);
plot(t,PAM_Modulated,'linewidth',1.5);
title('4 PAM Modulated Signal')
ylim([-4 4])
linkaxes([ax1,ax2],'x');
grid on;


% Average Energy of 4-PAM - Evaluation
M=4;
fun = @(t) (3*cos(2*pi*Fcarrier*t/fsamp)).^2;
q1 = integral(fun,0,1/fsamp);

fun = @(t) (cos(2*pi*Fcarrier*t/fsamp)).^2;
q2 = integral(fun,0,1/fsamp);

Es_avg_4PAM_Evaluation=1/M*(2*q1+2*q2)


% Average Energy of 4-PAM - Formula
Es_avg_4PAM_Formula=Eg*(M^2-1)/3



disp('****************************************************************')






% (((((((((((((((((( 4-FSK )))))))))))))))))))))

FSK_Modulated=[];

t_samp=0:1/fsamp:2-1/fsamp;

for k=1:2:length(b)
    if (b(k)==0 && b(k+1)==0)
        S_M=cos(2*pi*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==1 && b(k+1)==0)
        S_M=cos(2*pi*2*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==0 && b(k+1)==1)
        S_M=cos(2*pi*3*Fcarrier*t_samp/fsamp);
    end
    if (b(k)==1 && b(k+1)==1)
        S_M=cos(2*pi*4*Fcarrier*t_samp/fsamp);
    end
    FSK_Modulated=[FSK_Modulated S_M];
end

figure(4)
ax1 =subplot(2,1,1);
plot(t,BPSK_Baseband,'linewidth',1.5);
ylim([-2 2])
xlabel('$t/T_s$','interpreter','latex');
ylabel('$s(t)$','interpreter','latex');
title('BPSK BaseBand Signal')
grid on;
ax2 =subplot(2,1,2);
plot(t,FSK_Modulated,'linewidth',1.5);
title('4 FSK Modulated Signal')
ylim([-2 2])
linkaxes([ax1,ax2],'x');
grid on;


% Average Energy of 4-FSK - Evaluation
M=4;
fun = @(t) (cos(2*pi*Fcarrier*t/fsamp)).^2;
q1 = integral(fun,0,1/fsamp);
fun = @(t) (cos(2*pi*2*Fcarrier*t/fsamp)).^2;
q2 = integral(fun,0,1/(2*fsamp));
fun = @(t) (cos(2*pi*3*Fcarrier*t/fsamp)).^2;
q3 = integral(fun,0,1/(3*fsamp));
fun = @(t) (cos(2*pi*4*Fcarrier*t/fsamp)).^2;
q4 = integral(fun,0,1/(4*fsamp));

Es_avg_4FSK_Evaluation=(1/M)*(q1+q2+q3+q4)


% Average Energy of 4-FSK - Formula
Es_avg_4FSK_Formula=Eg/2






disp('****************************************************************')






% (((((((((((((((((( 16-QAM )))))))))))))))))))))

% Amplitude and Phase of 16-QAM
index=0;
QAM16_MAT=[];
for k=3:-2:-3
    for l=-3:2:3
        QAM16_MAT=[QAM16_MAT; index abs(l+k*i) angle(l+k*i)];
        index=index+1;
    end
end



QAM16_Modulated=[];

t_samp=0:1/fsamp:4-1/fsamp;

for k=1:4:length(b)
    if (b(k)==0 && b(k+1)==0 && b(k+2)==0 && b(k+3)==0)
        symbol_index=1;
    end
    if (b(k)==1 && b(k+1)==0 && b(k+2)==0 && b(k+3)==0)
        symbol_index=2;
    end
    if (b(k)==0 && b(k+1)==1 && b(k+2)==0 && b(k+3)==0)
        symbol_index=3;
    end
    if (b(k)==1 && b(k+1)==1 && b(k+2)==0 && b(k+3)==0)
        symbol_index=4;
    end
    if (b(k)==0 && b(k+1)==0 && b(k+2)==1 && b(k+3)==0)
        symbol_index=5;
    end
    if (b(k)==1 && b(k+1)==0 && b(k+2)==1 && b(k+3)==0)
        symbol_index=6;
    end
    if (b(k)==0 && b(k+1)==1 && b(k+2)==1 && b(k+3)==0)
        symbol_index=7;
    end
    if (b(k)==1 && b(k+1)==1 && b(k+2)==1 && b(k+3)==0)
        symbol_index=8;
    end  
    if (b(k)==0 && b(k+1)==0 && b(k+2)==0 && b(k+3)==1)
        symbol_index=9;
    end
    if (b(k)==1 && b(k+1)==0 && b(k+2)==0 && b(k+3)==1)
        symbol_index=10;
    end
    if (b(k)==0 && b(k+1)==1 && b(k+2)==0 && b(k+3)==1)
        symbol_index=11;
    end
    if (b(k)==1 && b(k+1)==1 && b(k+2)==0 && b(k+3)==1)
        symbol_index=12;
    end
    if (b(k)==0 && b(k+1)==0 && b(k+2)==1 && b(k+3)==1)
        symbol_index=13;
    end
    if (b(k)==1 && b(k+1)==0 && b(k+2)==1 && b(k+3)==1)
        symbol_index=14;
    end
    if (b(k)==0 && b(k+1)==1 && b(k+2)==1 && b(k+3)==1)
        symbol_index=15;
    end
    if (b(k)==1 && b(k+1)==1 && b(k+2)==1 && b(k+3)==1)
        symbol_index=16;
    end       
    S_M=QAM16_MAT(symbol_index,2)*cos(2*pi*Fcarrier*t_samp/fsamp+QAM16_MAT(symbol_index,3));
    QAM16_Modulated=[QAM16_Modulated S_M];
end

figure(5)
ax1 =subplot(2,1,1);
plot(t,BPSK_Baseband,'linewidth',1.5);
ylim([-2 2])
xlabel('$t/T_s$','interpreter','latex');
ylabel('$s(t)$','interpreter','latex');
title('BPSK BaseBand Signal')
grid on;
ax2 =subplot(2,1,2);
plot(t,QAM16_Modulated,'linewidth',1.5);
title('16 QAM Modulated Signal')
ylim([-6 6])
linkaxes([ax1,ax2],'x');
grid on;



% Average Energy of 16 QAM - Evaluation
M=16;

Q_sum=0;
for k=1:M
    fun = @(t) (QAM16_MAT(k,2)*cos(2*pi*Fcarrier*t/fsamp+QAM16_MAT(k,3))).^2;
    q = integral(fun,0,1/fsamp);
    Q_sum=Q_sum+q;
end

Es_avg_16QAM_Evaluation=1/M*(Q_sum)


% Average Energy of 16 QAM - Formula
Es_avg_16QAM_Formula=2*(M-1)/3*(Eg/2)






NFFT=2048;
s=(BPSK_Baseband+1)/2;
[R,f]=periodogram(s,[],NFFT,fsamp);
figure(6)
plot(f,10*log10(R))
grid on
title('Periodigram of  s(t)')
xlabel('Frequency (Hz)')
ylim([-50 20])

[R,f]=periodogram(BPSK_Modulated,[],NFFT,fsamp);
figure(7)
plot(f,10*log10(R))
grid on
title('Periodigram of  BPSK Modulated')
xlabel('Frequency (Hz)')
ylim([-50 20])

[R,f]=periodogram(QPSK_Modulated,[],NFFT,fsamp);
figure(8)
plot(f,10*log10(R))
grid on
title('Periodigram of  QPSK Modulated')
xlabel('Frequency (Hz)')
ylim([-50 20])


[R,f]=periodogram(PAM_Modulated,[],NFFT,fsamp);
figure(9)
plot(f,10*log10(R))
grid on
title('Periodigram of  PAM Modulated')
xlabel('Frequency (Hz)')
ylim([-50 20])


[R,f]=periodogram(FSK_Modulated,[],NFFT,fsamp);
figure(10)
plot(f,10*log10(R))
grid on
title('Periodigram of  FSK Modulated')
xlabel('Frequency (Hz)')
ylim([-70 20])



[R,f]=periodogram(QAM16_Modulated,[],NFFT,fsamp);
figure(11)
plot(f,10*log10(R))
grid on
title('Periodigram of  QAM16 Modulated')
xlabel('Frequency (Hz)')
ylim([-50 20])