clc
clear
close all


numSymb=100;
fsamp=100; % number samples per Ts
time=0:1/fsamp:numSymb-1/fsamp; % normalize by sympol duration Ts

fc=4*fsamp;

EbN0dB=2;
EbN0=10^(EbN0dB/10);


% random binary data
b=randi(2,1,numSymb)-1;
% A=2*B-1;

Data_Baseband=[];
for k=1:length(b)
    if (b(k)==0)
        S=zeros(1,fsamp);
    else
        S=ones(1,fsamp);
    end
    Data_Baseband=[Data_Baseband S];
end



% (((((((((((((((((( 4-QAM Modulated Signal)))))))))))))))))))))

% Amplitude and Phase of 4-QAM
index=0;
QAM4_MAT=[];
for k=1:-2:-1
    for l=1:-2:-1
        QAM4_MAT=[QAM4_MAT; index abs(l+k*i) angle(l+k*i)];
        index=index+1;
    end
end



QAM4_Modulated=[];
Symbols=[];
t_samp=0:1/fsamp:2-1/fsamp;

for k=1:2:length(b)
    if (b(k)==0 && b(k+1)==0)
        symbol_index=1;
    end
    if (b(k)==1 && b(k+1)==0)
        symbol_index=2;
    end
    if (b(k)==0 && b(k+1)==1)
        symbol_index=3;
    end
    if (b(k)==1 && b(k+1)==1)
        symbol_index=4;
    end
    S_M=QAM4_MAT(symbol_index,2)*cos(2*pi*fc*t_samp/fsamp+QAM4_MAT(symbol_index,3));
    QAM4_Modulated=[QAM4_Modulated S_M];
    Symbols=[Symbols symbol_index];
end

figure(1)
ax1 =subplot(4,1,1);
plot(time,Data_Baseband,'linewidth',1.5);
title('BaseBand Data')
ylim([-2 2])
grid on;
ax2 =subplot(4,1,2);
x=2:2:numSymb;
vals=Symbols;
bar1 = bar(x,vals,0.4);
xtips2 = bar1(1).XEndPoints;
xlim([0 numSymb])
ytips2 = bar1(1).YEndPoints;
labels2 = string(bar1(1).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0,4])


ax3 =subplot(4,1,3);
plot(time,QAM4_Modulated,'linewidth',1.5);
grid on;
title('4 QAM Modulated Signal')

z=QAM4_Modulated;

%AWGN channel
Es=sum(z.*z)/(numSymb*fsamp);         % average symbol energy
k=1; % binary
sigma2=fsamp*Es/(2*k*EbN0);
n=sqrt(sigma2)*randn(1,length(time));
r=z+n;


ax4 =subplot(4,1,4);
plot(time,r,'linewidth',1.5);
title('4 QAM Noisy Modulated Signal')
linkaxes([ax1,ax2,ax3,ax4],'x');
grid on;


% (((((((((((((((((( 4-QAM Receiver)))))))))))))))))))))

t=0:1/fsamp:2-1/fsamp;
detected_symbols=[];
for k=1:length(b)/2

    r1=r(2*(k-1)*fsamp+1:2*k*fsamp);
    funcx=r1.*cos(2*pi*fc*t/fsamp);
    funcy=-r1.*sin(2*pi*fc*t/fsamp);

    x=trapz(funcx);
    y=trapz(funcy);

    AiBi=[1 1;
        -1 1
        1  -1
        -1 -1];
    sai=[];
    for n=1:4
        sai=[sai; AiBi(n,1)*x+AiBi(n,2)*y];
    end
    detected_symbols=[detected_symbols find(sai==max(sai))];
end

SemmbolErrorRate=length(find(detected_symbols-Symbols)~=0)/length(Symbols)


M=4;
c=M*(1-1/sqrt(M));
dminSqure=3*log2(M)/(M-1);

Ps_Formulla=c*QFunc(sqrt(dminSqure*EbN0))

function error=QFunc(a)

error=0.5-0.5*erf(a/sqrt(2));

end



