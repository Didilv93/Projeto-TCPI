s = tf('s');
load('params.mat');
load('eqs_ex2.mat');

%%
p = 10; % Horizonte de Predicao
m = 5;  % Horizonte de Controle
lambda = 1;
N = 2500; 

%% Fiding parameters
% For power transfer function
[g, t] = lsim(Gz_p, ones(N+1, 1), T*(1:N+1)');
g = [g(2:end); g(end)];
G = triu(toeplitz(g(1:p)))';
H = inv(G'*G+lambda*eye(p))*G';
K = H(1, :);
[num, dem] = tfdata(Gz_p);
B = cell2mat(num);
A = cell2mat(dem);

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
w = w(1:2:end); % Obtem w 2x mais rapido q g

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
for t=1:N-T
    % Calcula y_t    
    y_t = -A(2:end)*flip(y(end-na+2:end))+B*flip(u(end-nb+1:end));
    
    % Salva y_t
    y(t) = y_t;
    
    %Calcula f
    f=[];
    for k=1:p
        f(k, 1) = y_t;
        for i=1:p
            f(k, 1) = f(k, 1) + ( g(k+i)-g(i) )*du;
        end
    end    

    
    % Resovler com restricoes
    w_p = w(t+1:t+p);
    
    H = 2*(delta*G'*G+eye(length(G)));
    b_T = 2*delta*(f-w_p)'*G;
    f0 = (f-w_p)'*(f-w_p);
    solution = quadprog(H, b_T, [], [], [], [], -0.6, 0.4);
    u(t+1) = solution(1);

    
    % Limita
    if(u(t+1) > 0.4)
        u(t+1) = 0.4;
    elseif (u(t+1)< -0.6)
        u(t+1) = -0.6;
    end
end


y_dmc_restrito = y;

figure();
stairs(1:N-T, y);
hold on; 
plot(w_plot);
plot(y_plot);
