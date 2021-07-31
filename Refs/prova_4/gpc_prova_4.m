s = tf('s');
load('params.mat');
load('eqs_ex2.mat');

%%
p = 3; % Horizonte de Predicao
m = 5;  % Horizonte de Controle
lambda = 1;
N = 2500; 

%% Fiding parameters
% For power transfer function
[g, t] = lsim(Gz_p, ones(N+1, 1), T*(1:N+1)');
G = [0 0 0; 0.0034 0 0; 0.0068 0.0034 0];
H = inv(G'*G+eye(length(G)))*G';
K = H(1, :);
fj = [0.0034 -0.0068 -0.0102];
[num, dem] = tfdata(Gz_p);
A = cell2mat(num);
B = cell2mat(dem);


F = [-19966 0.9966; 2.9898 -1.9898; 3.9593 2.9796];

% For disturbance transfer function
[d, t] = lsim(Gz_delta, ones(N+1, 1), T*(1:N+1)');
D = triu(toeplitz(d(1:p)))';
[num, dem] = tfdata(Gz_delta);
Bp = cell2mat(num);
Ap = cell2mat(dem);

%%
na = length(A);
nb = length(B);

%% W
[w, td] = lsim(Gz_p, ones(2*(N+1), 1), T*(1:2*(N+1))');
w = w(1:2:end);

%% Plot y vs W
[y_plot, t] = lsim(Gz_p, ones(1500, 1), T*(1:1500)');
[w_plot, td] = lsim(Gz_p, ones(2*1500, 1), T*(1:2*1500)');
w_plot = w_plot(1:2:end);
plot(t, y_plot);
hold on;
plot(t, w_plot);
legend('y', 'w');
t_teste = t;
%%

u = zeros(5, 1);
y = zeros(5, 1);
du = 0 ;
y_t = 0;
y = zeros(5, 1);
u = zeros(5 , 1);
y_alt = [];
for t=5:N-T
    % Calcula y_t    
    y_t = -A(2:end)*flip(y(end-na+2:end))+B*flip(u(end-nb+1:end));

        % Salva y_t
    y(t) = y_t;
    
    %Calcula f
    f=[];
    for k=1:p
        for i=1:length(F)
            f(k, 1) = fj(k)*du + F(k, :)*flip(y(end-1:end));
        end
    end    

    % Atualiza 
    du = K*(w(t+1:t+p) - f);
    u(t+1) = u(t) + du;
end

figure();
plot(1:N-T, y);
hold on; 
plot(w_plot)
