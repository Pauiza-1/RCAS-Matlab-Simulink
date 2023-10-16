yalmip('clear')
clear; 
clc;
format compact
format short e
%% Physical parameters
mb = 300;    % kg
mw = 60;     % kg
bs = 1000;   % N/m/s
ks = 16000 ; % N/m
kt = 190000; % N/m

% State matrices
A = [ 0         0      1       0;
      0         0      0       1;
      -ks/mb -ks/mb -bs/mb bs/mb; 
      ks/mw (-ks-kt)/mw bs/mw -bs/mw];

B2 = [ 0; 0; 1/mb;-1/mw]; % input matrix 1.7e4/ms in the B matrix improves the result. Why???
B1 = [ 0; 0; kt/mb;kt/mw];% disturbance matrix
C1 = [1 0 0 0;0 0 1 0]; % y measurement output
C2 = [0 1 0 0;0 0 0 1]; % z performance error output 
D = [0;0];

%% Form State Model 
%
qcar = ss(A,B2,C1,D);
% qcar.StateName = {'body travel (m)';'body vel (m/s)';...
%           'wheel travel (m)';'wheel vel (m/s)'};
% qcar.InputName = {'r';'fs'};
% qcar.OutputName = {'xb';'sd';'ab'};

Co = rank(ctrb(A,B2)); unco_sys = length(A) - Co; % Checks controllability

disp('No. of uncontr state = '); disp (unco_sys);

Ob1 = obsv(A,C1); unob1 = length(A)-rank(Ob1);
disp('No. of unob state = '); disp (unob1);

%% define variables
n=4; m=1; q=2;
gamma = sdpvar(m);
X = sdpvar(n,n);
Y = sdpvar(n,n);

A_hat = sdpvar(n,n);
B_hat = sdpvar(n,q);
C_hat = sdpvar(m,n);
N = eye(n); %sdpvar(n,n);%
M = sdpvar(n,n);
I = eye(n);
I2 = eye(m);
zero = zeros(q,m);

B1_r = [B1 zeros(4,3)]; % B1 all rows in bracket expression
C1_r = [X*C1' zeros(4,3)]; % C1 all rows in bracket expression

P = [X I;I Y];

Theta = [A*X+X*A'+B2*C_hat+C_hat'*B2' A_hat+A B1 X*C1';
       (A_hat+A)' Y*A+A'*Y+B_hat*C1+C1'*B_hat' Y*B1 C2';
       B1' (Y*B1)'  -gamma*I2  zero';
       (X*C1')' C2 zero -eye(q,q)*gamma];


constrants = [X>=0,Y>=0,P>=0,Theta<=0, M>=0];
   
options= sdpsettings('solver','mosek','verbose',3);
diagnostics = optimize(constrants,gamma,[]);
disp(diagnostics.problem);
if diagnostics.problem == 0
disp('Solution is Feasible')
check(constrants)
value(gamma)
elseif diagnostics.problem == 1
disp('Infeasible')
else
disp('Something else happened')
end
%disp('Press any key to continue . . . . Check for feasibility before controller parameters')

M=(I-X*Y)*inv(N');

% %  Controllers parameters
   Cc = value(C_hat)*(inv(value(M)'));
% % 
   Bc = (inv(value(N)))* value(B_hat);
% %  
  Ac = (inv(value(N)))*(value(A_hat) - value(Y)*A*value(X)-value(N)*Bc*C1*value(X)-value(Y)*B1*Cc*value(M)')*(inv(value(M)'));
disp('The closed-loop stability is: ')
  
  eig([A B2*Cc;Bc*C1 Ac])

%% Road bumps and potholes input  
t = 0:0.0025:10;
roaddist = zeros(size(t));
roaddist(501:601) = 0.04*(1-cos(8*pi*t(501:601)));
roaddist(1501:1601) = -0.04*(1-cos(8*pi*t(1501:1601)));
%% Nominal H infinity Control
D11 = [0;0]; D12 = [0;0]; D21 = [0;0]; D22 = [0;0];
B = [B1 B2];
C = [C2; C1];
Dn = [D11 D12;D21 D22];
ncont = 1;
nmeas = 3; % vehicle body acceleration, wheel velocity and  
opts = hinfsynOptions('Method','LMI','Display','on');
Pn = ss(A,B,C,Dn);
[K,CL,gam,info] = hinfsyn(Pn,nmeas,ncont,opts);