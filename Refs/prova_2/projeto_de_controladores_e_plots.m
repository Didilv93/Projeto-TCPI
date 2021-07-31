load('params.mat');
s = tf('s');

%% H1
G = qin_max/A1/s;
C = G^(-1)* 1/(tempo_acomodacao_h1/20*s) ;
step(feedback(G, 1));
hold on;
step(feedback(C*G, 1));
title('Resposta ao degrau - Nível');
legend('malha fechada', 'malha fechada cotrolada');
%% T1
G_P = delta1*qin_max/(A1*h1*s + delta1*qin_max);
tau_g_p = 1.7529e+04;
tau_des = 0.8*tau_g_p;

C = (tau_g_p*s+1)/(tau_des*s);
step(feedback(C*G_P, 1));
hold on 
step(G_P);
title('Resposta ao Degrau - Temperatura');
legend('malha fechada controlada', 'malha aberta');

zero_sys = tf([17529 1], [1]);

%% H1 com atraso
atraso = tempo_acomodacao_h1/20;
G = exp(-atraso*s)*qin_max/A1/s;
C = 1/(qin_max/A1/s)* 1/((tempo_acomodacao_h1/20+atraso/2)*s) ;
step(feedback(G, 1));
hold on;
step(feedback(C*G, 1));
title('Resposta ao degrau - Nível com Atraso');
legend('malha fechada', 'malha fechada cotrolada');

