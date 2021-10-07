% Reprodu��o de :
%
%             ON�OFF Error Control Coding Scheme for Minimizing 
%             Tracking Error in Wireless Feedback Control Systems
%
%       Shingo Hattori, , Kentaro Kobayashi, Hiraku Okada, Masaaki Katayama 
%
%
%    IEEE TRANSACTIONS ON INDUSTRIAL INFORMATICS, VOL. 11, NO. 6, DECEMBER 2015

clear 
close all
clc

%%%%%%%%       Defini��o do sistema de controle

%Par�metros do p�ndulo - Table 1
mp = 0.016; %Mass of pendulum rod (kg)
lp = 0.2; %Length of pendulum rod (m)
ra = 0.2; %Length of horizontal arm (m)
Jb = 0.0048; %Moment of inertia of the arm (kg*m^2)
Rm = 8.3; % DC motor armature resistance (Ohms)
Kg = 0.023; % DC Motor Constant (N*m/A, V*s/rad)
Km = 7.5; %Gear Ratio (120:16)
g = 9.81;% Acceleration of Gravity  (m/s^2)

% Matriz de estados cont�nua - (39)
Ac = [0                             0   1   0; 
      0                             0   0   1; 
      2*g*((Jb + mp*ra^2)/(lp*Jb))  0   0   -(2*ra*Kg^2*Km^2)/(Rm*lp*Jb);
      -(mp*ra*g)/Jb                 0   0   ((Kg^2)*(Km^2))/Rm*Jb];

Bc = [  0;
        0;
        -(2*ra*Kg*Km)/(Rm*lp*Jb);
        (Kg*Km)/Rm*Jb];

n = length(Ac); %n - ordem do sistema
out = n;%(m do espaco de estados) - n�mero de sa�das (neste caso assumiu-se que a sa�da � o vetor de estado)
in = 1;%(r do espaco de estados)- n�mero de entradas
    
Cc = eye(n);
Dc = zeros(out,in);
    
sys = ss(Ac,Bc,Cc,Dc);

%Par�metros de Amostragem - Table II        
Ts = 0.25;  %Sampling Period (s)
Fs = 1/Ts;  %Sampling Rate (Hz)
vezes = 3000;

    EndTime = 300; %End time (s)
    nSimul = EndTime/Ts; %n�mero de pontos da Simula��o
    pad = 3; %numero de bits a ser adicionados em nSimul porque o controle preditivo calcula amostras al�m do fim da simula��o   

    
    %calcula o equivalente discreto do sistema
    sysD = c2d(sys,Ts); %a d�vida � se isso funciona ou precisamos usar (41)
    [A,B,C,D] = ssdata(sysD); %retorna as matrizes de estado

    %Define o controlador LQR - Final da p�gina 1416
    %reparar que a defini��o dos estados utilizadas no artigo est� diferente da
    %utilizada no apendice, entao as duas colunas do meio da matriz Q devem ser
    %invertidas
    Q = [10 0 0 0;
        0 0.5 0 0;
        0 0 2 0;
        0 0 0 0];
    R = 0.5;

    %Calcula o ganho �timo do controlador LQR
    Ndist=0;% %matriz de peso do disturbio no c�lculo da fun��o objetivo do controlador LQR (ver help da fun��o no matlab)
    [K,S,e]=dlqr(A,B,Q,R,Ndist);

    %calcula o setpoint do sistema - (35)
    phi = [];
    T = 10; %logo ap�s 35
    t = 0:Ts:(nSimul+pad-1)*Ts;
    for kn=1:nSimul/(T/Ts)
        phi = [phi pi*ones(1,((T/2)/Ts)) zeros(1,((T/2)/Ts))];
    end
    phi = [phi zeros(1,3)]; %adiciona mais dois instantes de simula��o porque o controlador preditivo est� sempre 2 instantes a frente

    %�ltimo par�grafo, primeira coluna da p�gina 1417
    MaxAngle = pi/2; 
    MinAngle = -pi/2;%�ngulo m�ximo do p�ndulo - se o �ngulo for maior que isso, o p�ndulo caiu
    %caso ideal
    %caso uncoded

    K = -K;                       
    %Existem 4 subsistemas poss�veis (se��o 2)
    %S1 - sistema est�vel:Controlador-Planta OK, Planta-Controlador OK

    S1 = [(A+B*K) zeros(n,n) zeros(n,in);eye(n) zeros(n,n) zeros(n,in);K zeros(in,n) zeros(in,in)];

    %S2 - sistema inst�vel:Controlador-Planta OK, Planta-Controlador em Outage
    S2 = [A B*K zeros(n,in); zeros(n,n) eye(n) zeros(n,in); zeros(in,n) K zeros(in,in)];

    %S3 - sistema inst�vel:Controlador-Planta em Outage, Planta-Controlador OK
    S3 = [A zeros(n,n) B;eye(n) zeros(n,n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];

    %S4 - sistema inst�vel:Controlador-Planta em Outage, Planta-Controlador em Outage
    S4 = [A zeros(n,n) B;zeros(n,n) eye(n) zeros(n,in); zeros(in,n) zeros(in,n) eye(in)];

    lambda_1 = max(abs(eig(S1)));
    lambda_2 = max([abs(eig(S2)); abs(eig(S3)); abs(eig(S4))]);

initial_rho = 0;
initial_phi = pi/2-pi/16;%pi/32;

p1c = [0 1e-3 5e-3 7e-3 1e-2 1.3e-2];
p1s = [0 1e-3 5e-3 7e-3 1e-2 1.3e-2];
q = zeros(length(p1c),nSimul);     
tplots = 1; %plot temporal series?

outages_pdf = zeros(length(p1c),nSimul+pad);

for s = 1:vezes
    outageff = rand(1,nSimul+pad); %sorteio de outage de feedforward
    outagefb = rand(1,nSimul+pad); %sorteio de outage de feedback
                 
    for o = 1:length(p1c)    
        outages_pdf(o,:) = outages_pdf(o,:) + (((outageff < p1c(o)) + (outagefb < p1s(o))) >= 1);
        
       %resetar vari�veis
        x = [initial_phi*ones(1,nSimul+pad);(initial_rho)*ones(1,nSimul+pad);zeros(1,nSimul+pad);zeros(1,nSimul+pad)]; %estado do sistema
        u = zeros(in, nSimul+pad); %esfor�os de controle calculados pelo LQR
        w = zeros(n, nSimul+pad); 
        down = ones(1,nSimul+pad);        
        down(1) = 0;
        Z = [x;w;u];
            for k=2:nSimul
                if down(k-1)==1
                    if(o==4)
                        o=o+0;
                    end
                    q(o,k) = q(o,k) + 1;  
                    down(k-1) = 0;
                   break
                else
                    if outageff(k) > p1c(o) && outagefb(k) > p1s(o) %u[k] is received at k
                          Z(:,k) = S1*Z(:,k-1);  
                    elseif outageff(k) < p1c(o) && outagefb(k) > p1s(o)
                          Z(:,k) = S3*Z(:,k-1); 
                    elseif outageff(k) > p1c(o) && outagefb(k) < p1s(o)
                          Z(:,k) = S4*Z(:,k-1); 
                    else
                          Z(:,k) = S2*Z(:,k-1); 
                    end 
                    if abs(Z(1,k)) > MaxAngle 
                        down(k) = 1;             
                        %�ltimo par�grafo, primeira coluna da p�gina 1417
                        %If the angle of the pendulum |?[k]| > pi/6, we assume that the
                        %pendulum has fallen. Once the pendulum falls, the simulation
                        %run is terminated.
                    else
                        down(k) = 0;            
                    end                    
                end
                %%%% Lado do controlador
                %calcula os esfor�os de controle para a transmissao corrente
            end
            if tplots == 1
                [sizeZm sizeZn] = size(Z);
                Z_coded((o-1)*sizeZm+1:(o)*sizeZm,(s-1)*sizeZn+1:(s)*sizeZn) = Z; %estado do sistema
                out_coded(o,(s-1)*sizeZn+1:(s)*sizeZn) = (((outageff <  p1c(o)) + (outagefb <  p1s(o))) > 0);
                running(o,(s-1)*sizeZn+1:(s)*sizeZn) = ~down; 
                outageff_coded(1,(s-1)*sizeZn+1:(s)*sizeZn) = outageff; 
                outagefb_coded(1,(s-1)*sizeZn+1:(s)*sizeZn) = outagefb; 
            end
    end       
end

%q(1:length(p1c),1:5) = zeros(length(p1c),5);
S = cumsum(q,2);

figure
plot(S'./vezes,'LineWidth',3) 
for o = 1:length(p1c)  
    l1{o} = sprintf('Outage = %0.1e', p1c(o));
end
legend(l1)
%axis([3 60 0 0.9e-2])

if tplots == 1
    figure
    hold on
    [sizeZm, sizeZn] = size(Z);
    t = 1:length(Z_coded);
    plot([0 length(Z_coded)],[MaxAngle MaxAngle],'b:','LineWidth',2) 
    plot([0 length(Z_coded)],[MinAngle MinAngle],'b:','LineWidth',2) 
    lgd{1} = sprintf('Limite Máximo');
    lgd{2} = sprintf('Limite Mínimo');
    lo = 2;
    for o = 1:length(p1c)  
        angle = Z_coded((o-1)*sizeZm+1,:);
        angle(~running(o,:)) = NaN;
        stairs(t, angle,'LineWidth',3)
        lgd{lo+o} = sprintf('Outage = %0.1e', p1c(o));
    end
    for o = 1:length(p1c)  
        angle = Z_coded((o-1)*sizeZm+1,:);
        outages = out_coded(o,:).*angle;
        outages(~running(o,:)) = NaN;
        outages(outages==0) = NaN;
        plot(t,outages,'rX','LineWidth',3)
    end
    lgd{lo+length(p1c)+1} = sprintf('Outage');
    legend(lgd)
    axis([1.2391e5 1.2391e5+length(down)-1 MinAngle*1.3 -1*MinAngle*1.3])
end


t1 = 1.2391e5:(1.2391e5+length(down)-1);
v1 = outageff_coded(1,(t1));
v2 = outagefb_coded(1,(t1));
v1(1,44)=0.009;

% filename = 'test.mat';
% save(filename)
%%
% salvo outage e modificando para dar outage seguidas
% modificado a outage de feedback and feedfoward

for s = 1:vezes
    outageff = v1; %sorteio de outage de feedforward
    outagefb = v2; %sorteio de outage de feedback
                 
    for o = 1:length(p1c)    
        outages_pdf(o,:) = outages_pdf(o,:) + (((outageff < p1c(o)) + (outagefb < p1s(o))) >= 1);
        
       %resetar vari�veis
        x = [initial_phi*ones(1,nSimul+pad);(initial_rho)*ones(1,nSimul+pad);zeros(1,nSimul+pad);zeros(1,nSimul+pad)]; %estado do sistema
        u = zeros(in, nSimul+pad); %esfor�os de controle calculados pelo LQR
        w = zeros(n, nSimul+pad); 
        down = ones(1,nSimul+pad);        
        down(1) = 0;
        Z = [x;w;u];
            for k=2:nSimul
                if down(k-1)==1
                    if(o==4)
                        o=o+0;
                    end
                    q(o,k) = q(o,k) + 1;  
                    down(k-1) = 0;
                   break
                else
                    if outageff(k) > p1c(o) && outagefb(k) > p1s(o) %u[k] is received at k
                          Z(:,k) = S1*Z(:,k-1);  
                    elseif outageff(k) < p1c(o) && outagefb(k) > p1s(o)
                          Z(:,k) = S3*Z(:,k-1); 
                    elseif outageff(k) > p1c(o) && outagefb(k) < p1s(o)
                          Z(:,k) = S4*Z(:,k-1); 
                    else
                          Z(:,k) = S2*Z(:,k-1); 
                    end 
                    if abs(Z(1,k)) > MaxAngle 
                        down(k) = 1;             
                        %�ltimo par�grafo, primeira coluna da p�gina 1417
                        %If the angle of the pendulum |?[k]| > pi/6, we assume that the
                        %pendulum has fallen. Once the pendulum falls, the simulation
                        %run is terminated.
                    else
                        down(k) = 0;            
                    end                    
                end
                %%%% Lado do controlador

            end
            if tplots == 1
                [sizeZm sizeZn] = size(Z);
                Z_coded((o-1)*sizeZm+1:(o)*sizeZm,(s-1)*sizeZn+1:(s)*sizeZn) = Z; %estado do sistema
                out_coded(o,(s-1)*sizeZn+1:(s)*sizeZn) = (((outageff <  p1c(o)) + (outagefb <  p1s(o))) > 0);
                running(o,(s-1)*sizeZn+1:(s)*sizeZn) = ~down; 
            end
    end       
end


t2 = 1:80;
if tplots == 1
    figure
    hold on
    [sizeZm, sizeZn] = size(Z);
    lo = 2;
    for o = 2:length(p1c)  
        ax(:,o-1) = subplot((length(p1c)-1),1,o-1);
        
        angle_ideal = Z_coded(1,:);
        angle_ideal(~running(1,:)) = NaN;
        ang_ideal = angle_ideal(:,(t2));
        
        angle = Z_coded((o-1)*sizeZm+1,:);
        angle(~running(o,:)) = NaN;
        ang_end = angle(:,(t2));
        periodo = 1:length(ang_end); 
        
        outages = out_coded(o,:).*angle;
        outages(~running(o,:)) = NaN;
        outages(outages==0) = NaN;
        
        stairs(periodo,ang_ideal,'LineWidth',2)
        hold on
        stairs(periodo,ang_end,':','LineWidth',2)
        lgd{1} = sprintf('Sem Outage');
        lgd{2} = sprintf('Outage = %0.1e', p1c(o));
        lgd{3} = sprintf('Outage');
        plot(periodo,outages(:,(t2)),'kX','LineWidth',3)
%         plot(periodo,ones(size(periodo))*(17*pi/36),'-.','LineWidth',2) 
%         lgd{4} = sprintf('Initial $\phi$');
        grid on
        hold off
    	leg = legend(lgd(:,1:3),'Location','best')
        set(leg,'Interpreter','latex');
        set(leg,'FontSize',10);
        %axis([0 80 -pi/2 pi/2])
        %xticks([-pi/2 0 pi/2])
        %yticklabels({'-\pi/2','0','\pi/2'})
    end
     
end

    [ax1,h1]=suplabel('Tempo [s]');
    [ax2,h2]=suplabel('\^Angulo do Bra\,co ($\phi$)','y');
    set(h1,'Interpreter','latex');
    set(h1,'FontSize',17);   
    set(h2,'Interpreter','latex');
    set(h2,'FontSize',17); 

