clear
close all
clc
%Configurar matlabpool para usar 4 threads default para o "parfor" do
%mapping (executar uma vez no Matlab e nunca mais):
% myCluster = parcluster('local');
% myCluster.NumWorkers = 6;
% saveProfile(myCluster);

%ON–OFF Error Control Coding Scheme for Minimizing 
%Tracking Error in Wireless Feedback Control Systems
%Parâmetros do pêndulo - Table 1
mp = 0.016; %Mass of pendulum rod (kg)
lp = 0.2; %Length of pendulum rod (m)
ra = 0.2; %Length of horizontal arm (m)
Jb = 0.0048; %Moment of inertia of the arm (kg*m^2)
Rm = 8.3; % DC motor armature resistance (Ohms)
Kg = 0.023; % DC Motor Constant (N*m/A, V*s/rad)
Km = 7.5; %Gear Ratio (120:16)
g = 9.81;% Acceleration of Gravity  (m/s^2)

%Matriz de estados contínua - (39)
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

%Começo do Programa
n = length(Ac);                                 %ordem do sistema
out = n;                                        %(m do espaco de estados) - número de saídas 
in = 1;                                         %(r do espaco de estados) - número de entradas
Cc = eye(n);
Dc = zeros(out,in);
%modelo de espaço de estados
sys = ss(Ac,Bc,Cc,Dc);                           
Ts = 0.25;                     %tempo de amostragem
Fs = 1/Ts;                     %taxa de amostragem
EndTime = 10;                  %End time (s) arrumar depois 
nSimul = EndTime/Ts;           %número de pontos da Simulação 
t = 0:Ts:(nSimul)*Ts;
%discretizar o sistema continuo para obtenção das matrizes discretas
sysD = c2d(sys,Ts);            %discretização do sistema
[A,B,C,D] = ssdata(sysD);      %retorna as matrizes de estado discreto
%Controlador ótimo(LQR) com o sistema discretizado
Q = [10 0 0 0;
    0 0.5 0 0;
    0 0 2 0;
    0 0 0 0];
R = 0.5;
N = 0;        
%  Q = eye(n);
%  R = 10;
%  N = 0;    
%matriz de peso do disturbio no cálculo da função objetivo do controlador LQR
[K,S,e]= dlqr(A,B,Q,R,N);      %quais os parâmetros utilizados para obter o valor de K do artigo?

%Lembrando que fizemos essa mudança de sinal do vetor pois no artigo a
%notação usada pelo autor é de u(k) = K*x(k), e a notação usada pelo MATLAB
%é u(k) = - K*x(k). Sendo assim é necessário fazer a mudança de sinal no vetor.
K = -K;                       
%K = [0.0399 0.0217 -0.8172] valor utilizado no artigo
%Existem 4 subsistemas possíveis (seção 2)
%S1 - sistema estável:Controlador-Planta OK, Planta-Controlador OK
%S1 = [(A+B*K);eye(3);K];
%S1 = [S1 zeros(7,4)];
S1 = [(A+B*K) zeros(n,n) zeros(n,in);eye(n) zeros(n,n) zeros(n,in);K zeros(in,n) zeros(in,in)];

%S2 - sistema instável:Controlador-Planta OK, Planta-Controlador em Outage
%S2 = [A;zeros(4,3)];
%S2 = [S2 [B*K; eye(3);K] zeros(7,1)];
S2 = [A B*K zeros(n,in); zeros(n,n) eye(n) zeros(n,in); zeros(in,n) K zeros(in,in)];

%S3 - sistema instável:Controlador-Planta em Outage, Planta-Controlador OK
%S3 = [A;eye(3);zeros(1,3)];
%S3 = [S3 zeros(7,3) [B;zeros(3,1); eye(1)]];
S3 = [A zeros(n,n) B;eye(n) zeros(n,n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];
%S4 - sistema instável:Controlador-Planta em Outage, Planta-Controlador em Outage
%S4 = [A;zeros(4,3)];
%S4 = [S4 [zeros(2,3);eye(3);zeros(2,3)] [B;zeros(3,1);eye(1)]];
S4 = [A zeros(n,n) B;zeros(n,n) eye(n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];

% lambda_1 se refere à matriz estável (S1) - máximo autovalor
% lambda_2 se refere às matrizes instáveis (S2, S3 e S4)
lambda_1 = max(abs(eig(S1)));
lambda_2 = max([abs(eig(S2)); abs(eig(S3)); abs(eig(S4))]);

k = 100000;
% Cálculo do h1,h2,h3,h4
% O Lemma 1 estabelece que no caso do sistema estável (S1), existe um
% escalar h_1, válido para qualquer k>=1, que atende:
% ||S_1^k|| =< h_1 lambda_1^k   ou   h_1 >= ||S_1^k||/lambda_1^k
% No caso dos sistemas instáveis (i = {2,3,4}):
% ||S_i^k|| =< h_i lambda_2^k   ou   h_i >= ||S_i^k||/lambda_2^k
% Portanto, criamos vetores com as normas para vários valores de k
% A normalização mais adequada e mais próxima aos valores do artigo foi a
% normalização do tipo 1. O valor do artigo utilizado para o maior h é
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

%Sabendo que r1 = taxa de sucesso; pode-se associar este mesmo parï¿½metro a
%1 - taxa de falha, faz-se a taxa de falha o valor da probabilidade de
%outage
%Usando o artigo da Goldsmith para pegar os valores que serï¿½o apresentados 
%abaixo, "Energy-Constrained Modulation Optimization".
%Table I - MQAM Parameters
%Fazendo a potï¿½ncia de transmissï¿½o variï¿½vel
PtdB = -99:100; 
Rb = 10e6;
Pt = 10.^(PtdB./10);
Ptx = 97.9e-3; %100mW
Prx = 112e-3; %100mW
NodB = -174-30; %dB
No = 10^(NodB./10); %No linear
B = 10e6; %B=10kHz
fc = 2.5e9; %fc=2.5GHz
Gl = 5; %Gl = 30dB
Gl_n = 10^(Gl./10);
Ml = 20; %Ml = 20dB
Ml_n = 10^(Ml./10);
n = 0.35;

%--------------------------------------------------------
%Calculo da potï¿½ncia recebida na referencia de 1m
Pr = (Pt.*Gl_n*(3e8/fc)^2)/((4*pi)^2);
%distï¿½ncia variando atï¿½ 100m, log distancia exp 4
laco = 1:100;
begin = lambda_1;
inicio = (begin + 0.01);
passo = ((1-inicio)/length(PtdB));
lambda = inicio:passo:1-passo;

for d = 1:length(laco)
        Prd = (Pr.*(1./laco(d)).^4)/Ml_n;
        %EbNo mï¿½nima
        EbNo = Prd./(No*Rb);
        x  = ((2^(Rb/B)) - 1)./EbNo ;
        %Nï¿½mero de antenas
        M2 = 1;     
        %Fazer com que a probabilidade de outage mï¿½xima seja = 10^-3
        %Colocar um if para se for maior fazer a probabilidade de outage ser igual
        %a 10^-3
        for y = 1:length(lambda)
            Pout(y,:) = gammainc(x,M2);
            Ptot(y,:) = (Pt./n) + Ptx + Prx;
        end
        
        pout_sc = Pout;
        pout_ca = Pout;
        r1 = (1-pout_sc).*(1-pout_ca);
        r2 = (1-pout_sc).*pout_ca;
        r3 = pout_sc.*(1-pout_ca);
        r4 = pout_sc.*pout_ca;
        r = (r2 + r3 + 2*r4)./2;
        %Lembrando temos k = sum(ni), e as taxas de eventos descritas para cada
        %subsistema pode ser definida como ri = ni/k
        %n = quantidade de pacotes que foram entregues por esse sistema.
        n1 = (r1).*k;
        n2  = (r2).*k;
        n3 = (r3).*k;
        n4 = (r4).*k;

        N_sigma = 2.*(n2+n3+n4);
        tau_bar = k./N_sigma;
        
        for y = 1:length(lambda) 
           lambda_ast = exp((log(lambda_2) - log(lambda_1)).*r) + log(lambda_1);
           for x = 1:length(Pt) 
                    lambda_asterisco(y,x) = min((lambda(y)-passo), max((lambda_1+passo), lambda_ast(y,x)));

                    %taxa de perda de pacote mï¿½xima(outage)
                    r_max(y,x) = (log(lambda_asterisco(y,x)) - (log(lambda_1)))./((log(lambda_2))-(log(lambda_1)));
                    tau_a_asterisco(y,x) = (log(h))./((log(lambda(y))) - (log(lambda_asterisco(y,x))));
                
                    maior = tau_bar.*Ts;
                    menor(y,x) = tau_a_asterisco(y,x).*Ts;
                    ma =  maior(y,x);
                    me = menor(y,x);
                    if ma >= me
                    else
                        Ptot(y,x) = inf;
                    end  
           end
            %argmax(eficiência energética) = argmax(phi/potência total)
            [min_Ptot(d,y), pos_efi]= min(Ptot(y,:),[],2);
            max_Pout(d,y) = Pout(d,pos_efi);
        end 
            %p_max é a potência de transmissão otimizada            
            [opt_Ptot(d,:), pos_efi_opt] = min(Ptot(y,:),[],2);
end  
            opt_Pt = (opt_Ptot -(Ptx + Prx))*n;
            
            lambda_certo = 0.999999;
            K1 = (1/4)*(4 - (log(lambda_1)/log(h)) + (log(lambda_2)/log(h)) - sqrt(((-8*log(lambda./lambda_1)/log(h)) + (-4 + (log(lambda_1)/log(h)) - (log(lambda_2)/log(h)))^2)));
            for e = 1:length(laco)
            Pt_1(e,:) = -((((2^(Rb/B)) - 1)*No*B*(laco(e)^4)*Ml_n*((4*pi)^2)./(log(1 - K1)*Gl_n*(3e8/fc)^2)));
            opt_ptot(e,:)=(Pt_1(e,:)./n) + Ptx + Prx;
            end
            
            
cor1 = [0, 0.4470, 0.7410];
cor2 = [0.8500, 0.3250, 0.0980];
cor3 = [0.9290, 0.6940, 0.1250];
cor4 = [0.4940, 0.1840, 0.5560];
cor5 = [0.4660, 0.6740, 0.1880];
cor6 = [0 0 0];
cor7 = [0.3010, 0.7450, 0.9330];

figure
semilogy(lambda,((opt_ptot(20,:)/Rb)),'color',cor7,'LineWidth',2)
hold all
semilogy(lambda(3:15:end),((min_Ptot(20,(3:15:end))/Rb)),'o','MarkerEdge',cor1,'MarkerFace',cor1,'MarkerSize',8,'LineWidth',2)
semilogy(lambda,((opt_ptot(30,:)/Rb)),'color',cor7,'LineWidth',2)
semilogy(lambda(3:15:end),((min_Ptot(30,(3:15:end))/Rb)),'s','MarkerEdge',cor2,'MarkerFace',cor2,'MarkerSize',8,'LineWidth',2)
semilogy(lambda,((opt_ptot(50,:)/Rb)),'color',cor7,'LineWidth',2)
semilogy(lambda(3:15:end),((min_Ptot(50,(3:15:end))/Rb)),'d','MarkerEdge',cor3,'MarkerFace',cor3,'MarkerSize',8,'LineWidth',2)
semilogy(lambda,((opt_ptot(70,:)/Rb)),'color',cor7,'LineWidth',2)
semilogy(lambda(3:15:end),((min_Ptot(70,(3:15:end))/Rb)),'*','MarkerEdge',cor4,'MarkerFace',cor4,'MarkerSize',8,'LineWidth',2)
semilogy(lambda,((opt_ptot(90,:)/Rb)),'color',cor7,'LineWidth',2)
semilogy(lambda(3:15:end),((min_Ptot(90,(3:15:end))/Rb)),'p','MarkerEdge',cor5,'MarkerFace',cor5,'MarkerSize',8,'LineWidth',2)


w1= plot(50,50,'color',cor7,'LineWidth',2);
w2= plot(50,50,'o','MarkerEdge',cor1,'MarkerFace',cor1,'MarkerSize',8,'LineWidth',2);
w3= plot(50,50,'color',cor7,'LineWidth',2);
w4= plot(50,50,'s','MarkerEdge',cor2,'MarkerFace',cor2,'MarkerSize',8,'LineWidth',2);
w5= plot(50,50,'color',cor7,'LineWidth',2);
w6= plot(50,50,'d','MarkerEdge',cor3,'MarkerFace',cor3,'MarkerSize',8,'LineWidth',2);
w7= plot(50,50,'color',cor7,'LineWidth',2);
w8= plot(50,50,'*','MarkerEdge',cor4,'MarkerFace',cor4,'MarkerSize',8,'LineWidth',2);
w9= plot(50,50,'color',cor7,'LineWidth',2);
w10= plot(50,50,'p','MarkerEdge',cor5,'MarkerFace',cor5,'MarkerSize',8,'LineWidth',2);

grid on

t1 = ((min_Ptot(95,1)/Rb));
t2 = ((min_Ptot(10,end)/Rb));
axis([0.893 1 8*10^-8 2*10^-3])

leg1 = xlabel('Taxa de Decaimento ($\lambda$)')
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15); 
leg2 = ylabel('Energia Consumida por Bit ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);  

leg3 = legend([w1 w2 w4 w6 w8 w10],{'Resultados Anal\`iticos','d = 20m','d = 30m','d = 50m','d = 70m','d = 90m'},'Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10);             
