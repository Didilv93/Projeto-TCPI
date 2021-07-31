s = tf('s');


%%

G_delta = qin_max*(Tf-T1)/(delta1*qin_max)/(A1*h1*s/(delta1*qin_max)+1);
G_tf = (1/rho/c)/(A1*h1*s+delta1*qin_max);
G_p = 1/(A1*h1*s/(delta1*qin_max)+1);

%%

T = 60;
tau = A1*h1/(delta1*qin_max);
%% 

Gz_delta = c2d(G_delta, T, 'zoh');
Gz_tf = c2d(G_tf, T, 'zoh');
Gz_p = c2d(G_p, T, 'zoh');

