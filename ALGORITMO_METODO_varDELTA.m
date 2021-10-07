%ALGORITMO FINAL PARA TODOS OS TIPOS DE SISTEMAS
clear
close all
clc
%Configurar matlabpool para usar 4 threads default para o "parfor" do
%mapping (executar uma vez no Matlab e nunca mais):
% myCluster = parcluster('local');
% myCluster.NumWorkers = 6;
% saveProfile(myCluster);

%ON�OFF Error Control Coding Scheme for Minimizing 
%Tracking Error in Wireless Feedback Control Systems
%Par�metros do p�ndulo - Table 1
mp = 0.016; %Mass of pendulum rod (kg)
lp = 0.2; %Length of pendulum rod (m)
ra = 0.2; %Length of horizontal arm (m)
Jb = 0.0048; %Moment of inertia of the arm (kg*m^2)
Rm = 8.3; % DC motor armature resistance (Ohms)
Kg = 0.023; % DC Motor Constant (N*m/A, V*s/rad)
Km = 7.5; %Gear Ratio (120:16)
g = 9.81;% Acceleration of Gravity  (m/s^2)

%Matriz de estados cont�nua - (39)
Ac = [0                             0   1   0; 
      0                             0   0   1; 
      2*g*((Jb + mp*(ra^2))/(lp*Jb))  0   0   -(2*ra*(Kg^2)*(Km^2))/(Rm*lp*Jb);
      -(mp*ra*g)/Jb                 0   0   ((Kg^2)*(Km^2))/Rm*Jb];

Bc = [  0;
        0;
        -(2*ra*Kg*Km)/(Rm*lp*Jb);
        (Kg*Km)/Rm*Jb];

%Stability of networked control systems with packet dropout: an average
%dwell time approach
%Y.Sun, S.Qin
% Ac = [-1 0 -0.5; 1 -0.5 0;0 0 0.5];
% Bc = [0;0;1];

%Come�o do Programa
n = length(Ac);                                 %ordem do sistema
out = n;                                        %(m do espaco de estados) - n�mero de sa�das 
in = 1;                                         %(r do espaco de estados) - n�mero de entradas
Cc = eye(n);
Dc = zeros(out,in);
%modelo de espa�o de estados
sys = ss(Ac,Bc,Cc,Dc);                           
Ts = 0.25;                     %tempo de amostragem
Fs = 1/Ts;                     %taxa de amostragem
EndTime = 10;                  %End time (s) arrumar depois 
nSimul = EndTime/Ts;           %n�mero de pontos da Simula��o 
t = 0:Ts:(nSimul)*Ts;
%discretizar o sistema continuo para obten��o das matrizes discretas
sysD = c2d(sys,Ts);            %discretiza��o do sistema
[A,B,C,D] = ssdata(sysD);      %retorna as matrizes de estado discreto
%Controlador �timo(LQR) com o sistema discretizado
Q = [10 0 0 0;
    0 0.5 0 0;
    0 0 2 0;
    0 0 0 0];
R = 0.5;
N = 0;        
%  Q = eye(n);
%  R = 10;
%  N = 0;    
%matriz de peso do disturbio no c�lculo da fun��o objetivo do controlador LQR
[K,S,e]= dlqr(A,B,Q,R,N);      %quais os par�metros utilizados para obter o valor de K do artigo?

%Lembrando que fizemos essa mudan�a de sinal do vetor pois no artigo a
%nota��o usada pelo autor � de u(k) = K*x(k), e a nota��o usada pelo MATLAB
%� u(k) = - K*x(k). Sendo assim � necess�rio fazer a mudan�a de sinal no vetor.
K = -K;                       
%K = [0.0399 0.0217 -0.8172] valor utilizado no artigo
%Existem 4 subsistemas poss�veis (se��o 2)
%S1 - sistema est�vel:Controlador-Planta OK, Planta-Controlador OK
%S1 = [(A+B*K);eye(3);K];
%S1 = [S1 zeros(7,4)];
S1 = [(A+B*K) zeros(n,n) zeros(n,in);eye(n) zeros(n,n) zeros(n,in);K zeros(in,n) zeros(in,in)];

%S2 - sistema inst�vel:Controlador-Planta OK, Planta-Controlador em Outage
%S2 = [A;zeros(4,3)];
%S2 = [S2 [B*K; eye(3);K] zeros(7,1)];
S2 = [A B*K zeros(n,in); zeros(n,n) eye(n) zeros(n,in); zeros(in,n) K zeros(in,in)];

%S3 - sistema inst�vel:Controlador-Planta em Outage, Planta-Controlador OK
%S3 = [A;eye(3);zeros(1,3)];
%S3 = [S3 zeros(7,3) [B;zeros(3,1); eye(1)]];
S3 = [A zeros(n,n) B;eye(n) zeros(n,n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];
%S4 - sistema inst�vel:Controlador-Planta em Outage, Planta-Controlador em Outage
%S4 = [A;zeros(4,3)];
%S4 = [S4 [zeros(2,3);eye(3);zeros(2,3)] [B;zeros(3,1);eye(1)]];
S4 = [A zeros(n,n) B;zeros(n,n) eye(n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];

% lambda_1 se refere � matriz est�vel (S1) - m�ximo autovalor
% lambda_2 se refere �s matrizes inst�veis (S2, S3 e S4)
lambda_1 = max(abs(eig(S1)));
lambda_2 = max([abs(eig(S2)); abs(eig(S3)); abs(eig(S4))]);

k = 100000;
% C�lculo do h1,h2,h3,h4
% O Lemma 1 estabelece que no caso do sistema est�vel (S1), existe um
% escalar h_1, v�lido para qualquer k>=1, que atende:
% ||S_1^k|| =< h_1 lambda_1^k   ou   h_1 >= ||S_1^k||/lambda_1^k
% No caso dos sistemas inst�veis (i = {2,3,4}):
% ||S_i^k|| =< h_i lambda_2^k   ou   h_i >= ||S_i^k||/lambda_2^k
% Portanto, criamos vetores com as normas para v�rios valores de k
% A normaliza��o mais adequada e mais pr�xima aos valores do artigo foi a
% normaliza��o do tipo 1. O valor do artigo utilizado para o maior h �
% 14.0020

h1 = ones(1,k);
h2 = ones(1,k);
h3 = ones(1,k);
h4 = ones(1,k);

 for i = 1:k
       h1(i) = (norm(S1^i,1)./lambda_1.^i);
    if h1(i)<100
        x1 = max(h1);
    else
        break;
    end     
 end
     for i = 1:k
       h2(i) = (norm(S2^i,1)./lambda_2.^i);
    if h2(i)<100
        x2 = max(h2);
    else
        break;
    end
     end
  
    for i = 1:k
       h3(i) = (norm(S3^i,1)./lambda_2.^i);
    if h3(i)<100
        x3 = max(h3);
    else
        break;
    end
    end
    for i = 1:k
       h4(i) = (norm(S4^i,1)./lambda_2.^i);
    if h4(i)<100
        x4 = max(h4);
    else
        break;
    end
    end 
  
    h = max([x1 x2 x3 x4]);
%Sabendo que r1 = taxa de sucesso; pode-se associar este mesmo par�metro a
%1 - taxa de falha, faz-se a taxa de falha o valor da probabilidade de
%outage
%Usando o artigo da Goldsmith para pegar os valores que ser�o apresentados 
%abaixo, "Energy-Constrained Modulation Optimization".
%Table I - MQAM Parameters
%Fazendo a pot�ncia de transmiss�o vari�vel
PtdB = -40:0.01:10; 
Pt = 10.^(PtdB./10);
Ptx = 97.9e-3; %100mW
Prx = 112e-3; %100mW
NodB = -174-30; %dB
No = 10^(NodB./10); %No linear
B = 10e3; %B=10kHz
fc = 2.5e9; %fc=2.5GHz
Gl = 5; %Gl = 30dB
Gl_n = 10^(Gl./10);
Ml = 20; %Ml = 20dB
Ml_n = 10^(Ml./10);
n = 0.35;
alpha = 4;
c = 3e8;

%--------------------------------------------------------
%-------------------OTIMIZA��O TELECOMUNICA��ES-------


%-------------------RAYLEIGH-------------------
%dist�ncia variando at� 100m, log distancia exp 4

delta = 0.1:0.1:2.2;
laco = [10 15 20 25 30];
lambda = 0.999;
epsilon = 0.0001;
O_STAR = 1 - (log(lambda_1/lambda_2)/(4*log(h))) -(1/4)*(sqrt(((log(lambda_1/lambda_2)/(log(h))) -4)^2 - (8*log(lambda/lambda_1)/(log(h)))));
P_MAX = 10.^(10./10);
c = 3e8;

for d = 1:length(laco)
        TT = -(laco(d)^(-alpha))*(16*B*(laco(d)^alpha)*fc^2*Ml_n*No*(pi^2) + c^2*Gl_n*n*log(1-O_STAR)*(Ptx + Prx))/...
            (16*B*exp(1)*fc^2*Ml_n*No*(pi^2));
        RB = (1 + lambertw(TT))/log(2);
        PT = - (((4*pi)^2)*(No*B)*Ml_n*(laco(d)^alpha)*(2^(RB) - 1))/(log(1-O_STAR)*Gl_n*(c/fc)^2);
        PT = min(PT,P_MAX);

        PT_END(d,:) = 10*log10(PT);
        RB_END(d,:) = RB;
        EB_END(d,:) = ((PT/n) + Ptx + Prx)/(RB);
        
        EB_END_1 = EB_END;
end
    
for j = 1:length(delta) 
    
    laco = delta(j)*laco;
    for d = 1:length(laco)
            TT = -(laco(d)^(-alpha))*(16*B*(laco(d)^alpha)*fc^2*Ml_n*No*(pi^2) + c^2*Gl_n*n*log(1-O_STAR)*(Ptx + Prx))/...
                (16*B*exp(1)*fc^2*Ml_n*No*(pi^2));
            RB = (1 + lambertw(TT))/log(2);
            PT = - (((4*pi)^2)*(No*B)*Ml_n*(laco(d)^alpha)*(2^(RB) - 1))/(log(1-O_STAR)*Gl_n*(c/fc)^2);
            PT = min(PT,P_MAX);

            PT_END(d,:) = 10*log10(PT);
            RB_END(d,:) = RB;
            EB_END(d,:) = ((PT/n) + Ptx + Prx)/(RB);
    end
    
   EB_END_2(:,j) = EB_END;
   sum_EB_END(:,j) = EB_END_1 + EB_END_2(:,j);
end
    

figure
plot(delta,(sum_EB_END(1,:)),'-d','LineWidth',2)
hold all
plot(delta,(sum_EB_END(2,:)),'--','LineWidth',2)
plot(delta,(sum_EB_END(3,:)),'LineWidth',2)
plot(delta,(sum_EB_END(4,:)),'-s','LineWidth',2)
plot(delta,(sum_EB_END(5,:)),'-.','LineWidth',2)
axis([0 2.3 0.02 0.1])
grid on
leg1 = xlabel('Dist\^ancia Relativa entre Enlaces de Comunica\,c\~ao ($\delta$)');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15);      
leg2 = ylabel('Energia Consumida por Bir ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);   
leg3 = legend('$d_{sc} = 10\mathrm{m}$','$d_{sc} = 15\mathrm{m}$','$d_{sc} = 20\mathrm{m}$','$d_{sc} = 25\mathrm{m}$','$d_{sc} = 30\mathrm{m}$','Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10); 