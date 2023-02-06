%%
clear all;close all;
%% Requirements
[~,pList] = matlab.codetools.requiredFilesAndProducts('main.m');
[{pList.Name}';{'Simscape Electrical'}]
%% System Parameters

% Physical constants
r = 0.07; %m;
S = 0.0152; %m

% Steady states
u1_ss = -4; %V (napetost vzbujanja prve crpalke v delovni tocki)
u2_ss = -3.5; %V (napetost vzbujanja druge crpalke v delovni tocki)
x1_ss = 0.7712; %V (napetost prvega senzorja v delovni tocki)
x2_ss = 3.2841; %V (napetost drugega senzorja v delovni tocki)
x3_ss = 4.8098; %V (napetost tretjega senzorja v delovni tocki)
% oziroma
h1_ss = 0.3030; % m (normalna visina vode v prvem shranjevalniku)
h2_ss = 0.2290; % m  (normalna visina vode v drugem shranjevalniku)
h3_ss = 0.1520; % m  (normalna visina vode v tretjem shranjevalniku)


% Constants
Ks1 = -27.6937; %V/m
Ks2 = -27.5827; %V/m
Ks3 = -27.5087; %V/m
Kos1 = 9.1500; %V
Kos2 = 9.6000; %V
Kos3 = 9.0000; %V
K32 = -0.1107*10e-3; % m/s
K31 = 0.1747*10e-3; % m/s
K30 = 0.0360*10e-3; % m/s


%% Constrains
global u_min u_max
u_min = -10;
u_max = 10;
%% Linearized model
A = [-0.0125,  0.0126,  0.0000;
      0.0125, -0.0246,  0.0121;
      0.0000,  0.0120, -0.0212];
B = [-0.0091;
      0.0000;
      0.0000];
C = [0.0, 0.0, 1.0];
D = [0.0];
%x_initial = [x1_ss;
%             x2_ss;
%             x3_ss];
global model;
model = ss(A,B,C,D);
x_initial = [0;
             0;
             0];

maxTp = max(1./-pole(model))

model_firstOrder = tf(-1/maxTp, [1, 1/maxTp])
%% Simulate parameters
tsim = 60000; %s

%PAZI: To je dejansko napetost na IZHODNEM senzorju. Torej referenca je
%napetost na izhodnem senzorju, ki meri visino.
ref_ss = x3_ss;
d_ref = 0.5; %V



t1 = 40000; % Adaptation start time
t2 = 50000; % Adaptation end time

u = ref_ss + aprbs(t1,3*maxTp,2*d_ref)-d_ref;
d = diff(u);
idx = find(d) + 1;
idx = [1;idx];

ref_in = {};
% ref_in.Vals = [ref_ss, ref_ss+d_ref, ref_ss-d_ref, ref_ss];
% ref_in.Times = [0, 4000, 7000, 10000];
ref_in.Vals = [];
ref_in.Times = [];
for ii = 1:length(idx) - 1
     ref_in.Vals = [ref_in.Vals, u(idx(ii))];
     ref_in.Times = [ref_in.Times, idx(ii)];
end

% References after learning phase

ref_in.Vals = [ref_in.Vals, 4.400, 5.000, 5.200, 4.800];
ref_in.Times = [ref_in.Times, 42000, 44000, 45000, 48000];





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    REGULATOR - MIT     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gamma_f = 0.05;
gamma_q = 0.06;
sigma = 0.0005; % Forgetting factor: sigma modification
am = 0.008;
bm = 0.008;

refModel = tf(am, [1, bm]);

sign_b = -1;

%% Simulate


out = sim("model_mit.slx");
%% Results for validation

t = out.ref.Time;

%Reference
r = out.ref.Data;

% Regulator output signal
u = out.u.Data;

% Output
y = out.y.Data;

% Reference model output
yr = out.yr.Data;


J = out.J.Data(t>=t1 & t <=t2);
J = J(end)-J(1);
J


param_f = out.f.Data(end)
param_q = out.q.Data(end)


%%
figure;

subplot(2,1,1);
hold on;
plot(t, r);
plot(t, yr);
plot(t, y);
%ylim([0,2]);
xline(t1,'--');
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Voltage [V]', 'Interpreter', 'latex');
legend('$r_1$', '$y_{r}$',' $y_{p}$', 'Interpreter', 'latex');

subplot(2,1,2);
hold on;
plot(t, u);
yline([u_min u_max],'--');
xline(t1,'--');
ylim([-11, 11]);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Voltage [V]', 'Interpreter', 'latex');
legend('$u$', 'Interpreter', 'latex');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   REGULATOR - Lyapun   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gamma_f = 0.05;
gamma_q = 0.06;
sigma = 0.00005; % Forgetting factor: sigma modification
am = 0.008;
bm = 0.008;

refModel = tf(am, [1, bm]);

sign_b = -1;

%% Simulate


out = sim("model_lyapun.slx");

%% Results for validation

t = out.ref.Time;

%Reference
r = out.ref.Data;

% Regulator output signal
u = out.u.Data;

% Output
y = out.y.Data;

% Reference model output
yr = out.yr.Data;


J = out.J.Data(t>=t1 & t <=t2);
J = J(end)-J(1);
J


param_f = out.f.Data(end)
param_q = out.q.Data(end)


%%
figure;

subplot(2,1,1);
hold on;
plot(t, r);
plot(t, yr);
plot(t, y);
%ylim([0,2]);
xline(t1,'--');
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Voltage [V]', 'Interpreter', 'latex');
legend('$r_1$', '$y_{r}$',' $y_{p}$', 'Interpreter', 'latex');

subplot(2,1,2);
hold on;
plot(t, u);
yline([u_min u_max],'--');
xline(t1,'--');
ylim([-11, 11]);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Voltage [V]', 'Interpreter', 'latex');
legend('$u$', 'Interpreter', 'latex');



