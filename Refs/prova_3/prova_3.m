load('params.mat');

s = tf('s');


%% Ex 1  

g11 = Ro2*qq_max/(Ro2*A2*s+1);
g12 = Ro2*qf_max/(Ro2*A2*s+1);
g21 = (T1 - T2)*qq_max/(A2*h2*s + qq_max*deltaq +qf_max*deltaf);
g22 =  (Tf - T2)*qf_max/(A2*h2*s + qq_max*deltaq +qf_max*deltaf);

G = [g11 g12; g21 g22 ];

 G0 = [0.3766 0.8532; 25.3395 -34.0277];
 RGA_G = G0.*inv(G0)';
 
 %% Ex 2
 
 C1 = (2821*s+1)/(0.3766*0.8*2815*s);
 C2 = -(2820*s+1)/(34.0206*0.8*2820*s);
 
 %% Ex 3
 P = 1/( (2821*s+1)*(2820*s+1) );
 N = [ 0.3766*(2820*s+1)   0.8532*(2820*s+1);
      25.3395*(2821*s+1) -34.0277*(2821*s+1)];
  
 D = inv(N);
 
 %% Ex 4 
 ts = 0.8 * 3 * 2821;
 wn = 4 / 0.68 /ts;
 g_des = wn^2/(s^2+2*0.68*wn*s+wn^2);
 
Cd1 = (2821*s+1)*(2820*s+1)*7.5e-7/(s*(s+1.2e-3));

%% Ex 5
C_til = D*(Cd1*eye(2));
step(feedback(C_til*G,eye(2)));
c11 = C_til(1,1);
c12 = C_til(1,2);
c21 = C_til(2,1);
c22 = C_til(2,2);