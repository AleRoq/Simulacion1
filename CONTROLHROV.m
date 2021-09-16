
clear all; close all;
A=[-2.9595    0.1104    0.2761         0         0         0
   -0.0905   -1.9715   -0.3893         0         0         0
   -0.1673   -0.2855   -0.5433         0         0         0
    0.9950         0         0         0         0   -0.0429
    0.0998    1.0000         0         0         0    0.4279
         0         0    1.0000         0         0         0];
 
B=[ 0.1874    0.1874
    0.0001   -0.0001
   -0.3174    0.3174
         0         0
         0         0
         0         0];

C=[ 1     0     0     0     0     0
     0     1     0     0     0     0
     0     0     1     0     0     0
     0     0     0     1     0     0
     0     0     0     0     1     0
     0     0     0     0     0     1];

D=[  0     0
     0     0
     0     0
     0     0
     0     0
     0     0];

% Reducing the model
% X=[u, v, r]', U=[tau_u, tau_r]', Y=[u r]'
ap=A(1:3,1:3);  
bp=B(1:3,:);
cp=[C(1,1:3);
    C(3,1:3)];
dp=[D(1,1:2);
    D(3,1:2)];
    
[evec,eval] = eig(ap)   % evec contains eigenvectors
                        % eval contains poles or eigenvalues

t = [0:0.1:12];
u = [0*t' 0*t'];         % Set input u to zero for all time in order to generate zero input response;
                         % i.e. response to an initial condition x_o.

%
% Excite F404 (slow) temperature mode.
% This mode 
%            is associated with a pole at s = - 0.4.
%            has a time constant of  2.5 sec and settles in about 12.5 sec.
%            is the slowest of the F404's three (3) modes.
%            is associated with T_45
% Comment: It takes a long time to change temperature!
%
y = lsim(ss(ap, bp, eye(3,3), 0*ones(3,2)), u, t, evec(:,1)); 
plot(t,y)
grid
title('F404 Slow Temperature Mode: x_o = [ 0 0 1 ]')
ylabel('States')
xlabel('Time (seconds)')


y = lsim(ss(ap, bp, eye(3,3), 0*ones(3,2)), u, t, evec(:,3)); 
plot(t,y)
grid
title('F404 Slow RPM Mode: x_o = [ -0.9638 -0.2245 0.1440 ]')
ylabel('States')
xlabel('Time (seconds)')

%z = tzero(ss(ap,bp,cp,dp))                % transmission zeros
[p,z] = pzmap(ss(ap,bp,cp,dp))                % transmission zeros
zdir = null([z*eye(3)- ap  -bp; cp dp])   % transmission zero directions

% Controllability 
%
cm = [bp ap*bp (ap^2)*bp]  % Controllability Matrix
rcm= rank(cm)              % Rank of Controllability Matrix


%***************************************************************************
%
% Observability
%
om = [cp; 
      cp*ap;
      cp*(ap^2) ]          % Observability Matrix
rom = rank(om)             % Rank of Observability Matrix

%Ploteo de valores singulares de la planta

w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Outputs: N_2, T_45; Inputs: W_f/110, A_8/22')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

%Step 
sysG=ss(ap,bp,cp,dp)
step(sysG)

% Augment Plant with Integrators at Plant Input and Plot Singular Values
%
[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ]

b = [ 0*ones(ns,nc)
      eye(nc)      ]

c = [ cp  0*ones(nc,nc) ]

d = 0*ones(nc,nc)
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Design Plant Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

sv = sigma(ss(a, l, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Filter Open Loop (G_{FOL}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


pnint = eye(nc)                                    % Process Noise Intensity Matrix
mu = 0.01;                                         % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                   % Small mu - expensive sensor   - large bandwidth
                                                   % Large mu - inexpensive sensor - small bandwidth
mnint = mu*eye(nc)                                 % Measurement Noise Intensity Matrix 
% sysKal=ss(a, [b l], c, [d 0*ones(nc,nc)]);
% [kest, h, sig]= kalman(sysKal,pnint, mnint);  % Compute Filter Gain Matrix h
%[sig, poles, g1, rr] = care(a',c',l*l', mnint);                          
[sig, poles, g1] = care(a',c',l*l', mnint);  

                        
% Alternate Method for Computing h
h = g1';
sv = sigma(ss(a, h, c, d),w);
tsv = 20*log10(sv);
semilogx(w, tsv)
%clear sv
title('Target Loop (G_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

tolpoles = eig(a)                           % Target Open Loop Poles
%targzeros = tzero(a,h,c,0*ones(nc,nc))      % Target Open Loop Zeros
[targpoles,targzeros] = pzmap(ss(a,h,c,0*ones(nc,nc)))      % Target Open Loop Zeros
tclpoles = eig(a-h*c)                       % Target Closed Loop Poles

sv = sigma(ss(a-h*c, h, -c, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

sv = sigma(ss(a-h*c, h, c, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, 20*log10(10./w))
%clear sv
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                          % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(nc)                                      % Control Weigthing Matrix
%[k, poles, g, rr] = care(a,b,q,r);                   % Compute Control Gain Matrix G
[k, poles, g] = care(a,b,q,r);                   % Compute Control Gain Matrix G

% Compensator Analysis
%
ak = [ a-b*g-h*c  0*ones(ns+nc,nc)
       g          0*ones(nc,nc) ]

bk = [ h
       0*ones(nc,nc) ]

ck = [0*ones(nc, ns+nc) eye(nc,nc) ]

%cpoles = eig(ak)                               % Compensator Poles
%czeros = tzero(a, h, g, 0*ones(nc,nc))         % Compensator Zeros
[cpoles, czeros] = pzmap(ss(a, h, g, 0*ones(nc,nc)))         % Compensator Zeros
%zerocheck = tzero(ak, bk, ck, 0*ones(nc,nc))   % Check Compensator Zeros
[polecheck, zerocheck] = pzmap(ss(ak, bk, ck, 0*ones(nc,nc)))   % Check Compensator Zeros

sv = sigma(ss(ak, bk, ck, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Open Loop Analysis
%
al = [ ap                     bp*ck
       0*ones(ns+nc+nc,ns)    ak    ]

bl = [ 0*ones(ns,nc)
       bk ]
    
cl = [ cp  0*ones(nc,ns+nc+nc) ]
    
%olpoles = eig(al)                          % Open Loop Poles
%olzeros = tzero(al,bl,cl,0*ones(nc,nc))    % Open Loop Zeros
[olpoles, olzeros] = pzmap(ss(al,bl,cl,0*ones(nc,nc)))    % Open Loop Zeros    
sv = sigma(ss(al, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, tsv)
%clear sv
title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

% Closed Loop Analysis
%
clpoles = eig(al-bl*cl)           % Closed Loop Poles
clpkf = eig(a - h*c)              % Closed Loop Poles Due to Kalman Filter
clpreg = eig(a - b*g)             % Closed Loop Poles Due to Regulator


sv = sigma(ss(al-bl*cl, bl, -cl, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')

sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')


% Step Response in Closed Loop 
%

[y,t] = step(ss(al-bl*cl, bl, cl, 0*eye(nc)));
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')














