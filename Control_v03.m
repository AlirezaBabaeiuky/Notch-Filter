% Divide into 4 sections: 
% 1- no ref. tracking, only stability analysis for both collocated and
% non-colloated transfers - design compensator and notch filters for
% loop-shaping. this was indeed sections 1 and 2 
% sections 3 and 4: we have ref. tracking; then only focus on collocated
% (you can still define noncollocated wuth distrurbance at m1 and trasnfer
% measurement/assesment and ctrl designer at m2) but not worth the time and
% energy for now. So 3 and 4 are: do ref. tracking + notch or ref. tracking
% + compensator. of course for ref. tracking we need controllers + step
% response ... along with stability in bode/nyquist/root locus ...
% KEEP IT SIMPLE 
clc; 
clear all;
close all; 
m1 = 118;
m2 = 72;
k1 = 25e4;
k2 = 55e4; 
wn1_calc = sqrt(k1/m1)
wn2_calc = sqrt(k2/m2)
zeta1 = 1e-3;% 5e-5;
zeta2 = 8e-4;%25e-6;
ccr1 = 2*sqrt(k1*m1);
ccr2 = 2*sqrt(k2*m2);
c1 = zeta1*ccr1;
c2 = zeta2*ccr2; 

wca1 = sqrt(k1/(m1+m2)) % as rigid body mode freq. 
wca2 = sqrt(k2/m2) % as anti-resonance freq.
wca3 = sqrt(k2/(m1*m2/(m1+m2))) % as the second pole resonane there is 180 degrees of phase shift 

A_ss1 = [0, 1, 0, 0;
-(k1+k2)/m1, -(c1+c2)/m1, k2/m1, c2/m1; 
0, 0, 0, 1;
k2/m2, c2/m2, -k2/m2, -c2/m2]
B_ss1 = [0;
1/m1;
0;
0];
C_ss1 = [1, 0, 0, 0];
C_ssnc1 = [0, 0, 1, 0]; % Non-collocated measurement 
D_ss1 = 0; 

ssm1 = ss(A_ss1, B_ss1, C_ss1, D_ss1) % mss stnads for state-space model; both disturbance and measurement at m1 - Collocated transfer
ssmnc1 = ss(A_ss1, B_ss1, C_ssnc1, D_ss1) % mss stnads for state-space model; disturbance is at m1 and measurement at m2 - Non-colloate transfer
tfm1 = tf(ssm1)
tfmnc1 = tf(ssmnc1)

figure(1) 
bode(ssm1);
grid on 
hold on 
bode(ssmnc1);


[respm1, woutm1] = freqresp(ssm1);% frequency response by manual calc using: freqresp built-in fcn 
figure(2)
loglog(woutm1, abs(squeeze(respm1(1, 1, :))))
grid on;
figure(3)
semilogx(woutm1, angle(squeeze(respm1(1, 1, :))))
grid on;

%Notch Filtering - which is also a type of loop-shaping 
% Notch is also called: Band-Stop - Customized
wcnotch1 = 35; % Location of Notch. Notch cutoff / corner freq.
zetanotch1 = 0.01; % Depth of Notch - bigger value ends up in blunt/dull deptth
anotch1 = 0.05; % Width of Notch - bigger value ends up in narrower filter  
notchtfp1 = tf([1, 2*zetanotch1*wcnotch1, wcnotch1^2], [0, 0, wcnotch1^2]);
notchtfp2 = tf([0, 0, anotch1*wcnotch1], [1, anotch1*wcnotch1]);
notchtfp3 = tf([0, 0, wcnotch1/anotch1], [1, wcnotch1/anotch1]);
notchtfres1 = notchtfp1*notchtfp2*notchtfp3;% series(notchtfp1, notchtfp2, notchtfp3); 
% second formula for Notch filter - customized formula  
zetanotchnu = 0.2; % first of all; the numerator zeta should always be smaller than the denominator zeta;
% the bigger the difference, the deepe the notch!
zetanotchnures1nc = 0.0005; % for noncollocated resonance 1- using the same notch for the 2nd resonance 
zetanotchde = 0.1;
zetanotchderes1nc = 0.2; % for the noncolloated 
notchtf2 = tf([1, 2*zetanotchnu*wcnotch1, wcnotch1^2], [1, 2*zetanotchde*wcnotch1, wcnotch1^2]);
notchtf2res1nc = tf([1, 2*zetanotchnures1nc*wcnotch1, wcnotch1^2], [1, 2*zetanotchderes1nc*wcnotch1, wcnotch1^2]);

% Notch is also called: Band-Stop - Customized for 2nd resonance 
wcnotch2 = 115; % Location of Notch. Notch cutoff / corner freq.
zetanotch2 = 1e-3; % Depth of Notch - bigger value ends up in blunt/dull deptth
anotch2 = 0.002; % Width of Notch - bigger value ends up in narrower filter  
notchtfp21 = tf([1, 2*zetanotch2*wcnotch2, wcnotch2^2], [0, 0, wcnotch2^2]);
notchtfp22 = tf([0, 0, anotch2*wcnotch2], [1, anotch2*wcnotch2]);
notchtfp23 = tf([0, 0, wcnotch2/anotch2], [1, wcnotch2/anotch2]);
notchtfres21 = notchtfp21*notchtfp22*notchtfp23;% series(notchtfp1, notchtfp2, notchtfp3); 

notchtfall1 = notchtfres1 *  notchtfres21
notchtfall2 = notchtf2res1nc * notchtfres21

figure(4)
bode(notchtfall1)
hold on 
bode(notchtfall2)
grid on
hold off
% %legend('Customized notches', 'Standard notch and Customized notch', 'FontSize', 14)

% tfm1notchesCOL = tfm1 * notchtfall1 % tf for model 1 with notches Collocated measurement 
tfm1notchesCOL = tfm1 * notchtfall2 % tf for model 1 with notches Collocated measurement 
tfm1notchesNC = tfmnc1 * notchtfall2 % tf for model 1 with notches Non-Collocated measurement

figure(5)
bode(tfm1)
hold on 
bode(tfmnc1)
hold on 
bode(tfm1notchesCOL)
grid on 
hold on 
bode(tfm1notchesNC)
hold off
% %legend({'Freq. Resp. Collocated without notches', ...
%     'Freq. Resp. Non-Collocated without notches', 'Freq. Resp. Collocated with Customized Notches', ...
%     'Freq. Resp. Non-Collocated with Customized Notches'}, 'FontSize', 14)
% some notes: tuning the notch is also hard 

%% Section / project 3 (with ref. tracking and only collocated transfers) 
% BW is terrible idk, so go for Control1.m and look at the new design
% section and that is the project in there with good BW. 
m1 = 32;
m2 = 5;
k1 = 1e6;
k2 = 1e4; 
wn1_calc = sqrt(k1/m1)
wn2_calc = sqrt(k2/m2)
zeta1 = 5e-3;% 5e-5;
zeta2 = 25e-4;%25e-6;
ccr1 = 2*sqrt(k1*m1);
ccr2 = 2*sqrt(k2*m2);
c1 = zeta1*ccr1;
c2 = zeta2*ccr2;
m3 = 2.5; % mass of the spring added to the k2 
k3 = 1e9; % this way we shift the next res. freq. to very high values making 
% the sys as a ONE Rigid Body Mass 
mtot = m1+m2+m3;
A_ss2 = [0, 1;
-(k1+k2+k2)/mtot, -(c1+c2)/mtot] 
B_ss2 = [0;
1/mtot];
C_ss2 = [1, 0];
C_ssnc = [1, 0]; % Non-collocated measurement 
D_ss2 = 0; 
wca3 = sqrt(k1/mtot) 

ssm2 = ss(A_ss2, B_ss2, C_ss2, D_ss2) % mss stnads for state-space model

figure(22)
step(ssm2)
title('Note: this is Open loop step wo any ctrl ...')
xlim([0, 10]);
grid on 
StepInfossm2 = stepinfo(ssm2)
% No need to use pidTuner . I will tune the pid with editor command 
tfm2 = tf(ssm2)
% Design a notch for suppressing the peak at res freq. otherwise it can
% make sys unstable 
wcnotchss2 = 161; % Location of Notch. Notch cutoff / corner freq.
zetanotchss2 = 0.005; % Depth of Notch - bigger value ends up in blunt deptth
anotchss2 = 0.01; % Width of Notch - bigger value ends up in narrower filter  
notchtfp1ss2 = tf([1, 2*zetanotchss2*wcnotchss2, wcnotchss2^2], [0, 0, wcnotchss2^2]);
notchtfp2ss2 = tf([0, 0, anotchss2*wcnotchss2], [1, anotchss2*wcnotchss2]);
notchtfp3ss2 = tf([0, 0, wcnotchss2/anotchss2], [1, wcnotchss2/anotchss2]);
notchtfss2 = notchtfp1ss2*notchtfp2ss2*notchtfp3ss2;% series(notchtfp1, notchtfp2, notchtfp3); 

tfmnotchss2 = tfm2* notchtfss2;

tfcpid2 = pid(1e4, 5e6, 1e2);
tfcpid2 = pid(1e6, 5e6, 1e4);
tfcpid2 = pid(1e7, 5e6, 0); % more than 1e7 due to sat. limits on the mot is not
% possible. and also i increased BW as much as possible by increasing the
% Kp to have good BW then design compensators to address the bad ts
oltfm2notchpid = tfm2* notchtfss2*tfcpid2; 

figure(23) 
bode(ssm2) 
grid on 
hold on 
bode(tfmnotchss2) 
hold on 
bode(oltfm2notchpid) 
hold off 
% %legend({'New design MSD', 'MSD with Notch Filter', 'MSD with Notch and PID CTRL'}, 'FontSize', 16)

[Gmarss2, Pmarss2, WcroGss2, WcroPss2] = margin(oltfm2notchpid)
% recap: do all stability analysis like on bode / nyquist / Root Locus on OL and time
% response with Mp,ess, ts ... on CL 
cltfm2notchpid = feedback(oltfm2notchpid, 1, -1);
cltfm2notch = feedback(tfmnotchss2, 1, -1);
cltfm2 = feedback(tfm2, 1, -1);
figure(24)
step(cltfm2)
xlim([0, 10])
grid on 
hold on 
step(cltfm2notch)
hold on 
step(cltfm2notchpid)
% %legend({'new design', 'new design with notch', 'new design with Notch and PID'}, 'FontSize', 14)

[ycltfm2, tcltfm2] = step(cltfm2, 10);
[ycltfm2notch, tcltfm2notch] = step(cltfm2notch, 10);
[ycltfm2notchpid, tcltfm2notchpid] = step(cltfm2notchpid, 10);
essm2Per = (ycltfm2(end) - 1)*100
essm2notchPer = (ycltfm2notch(end) -1)*100 
essm2notchpidPer = (ycltfm2notchpid(end) -1)*100

SIAll1 = stepinfo(cltfm2)
SIAll2 = stepinfo(cltfm2notch)
SIAll3 = stepinfo(cltfm2notchpid)   % step info for all second design 
% Kp alone is terrible in ref track with significant e_ss so add I t oget
% rid of the offset 
% you can sue the pidTuner app but be very careful and import only and only
% CLTF  tfm2* notchtfss2;
save('tfm2.mat', 'tfm2')
save('notchtfss2.mat', 'notchtfss2')
% Sensitivity
Senoltfm2 = 1 / (1+ssm2)
Senoltfm2notch = 1 / (1+tfmnotchss2)
Senoltfm2notchpid = 1 / (1+oltfm2notchpid)
figure(241)
bode(Senoltfm2);
grid on 
hold on 
bode(Senoltfm2notch);
hold on 
bode(Senoltfm2notchpid);
hold off
% %legend({'Sensitivity of new design', 'Sensitivity of  new design with notch', 'Sensitivity of  new design with Notch and PID'}, 'FontSize', 14)
figure(242)
bode(Senoltfm2);
grid on 
% controlSystemDesigner(oltfm2notchpid) 
% pidTuner(oltfm2notchpid) % useless!!! unless you check for Bode only BUT NOT the stpe time 
% Does not need any compensators ... well design a compensator when you do
% not have any reference tracking ... 
% but indeed worked with pid as well so we can now design a compensator or
% better do fine-tuning on the pid 
Senss2 = 1 / (1+oltfm2notchpid)
figure(25)
bode(Senss2)
title("Sensitivity of 2nd Design - Accepetable Value of Sen is: mMAX of 5 dB ")
grid on 
figure(26)
rlocus(oltfm2notchpid)
grid on 
figure(27)
nyquist(oltfm2notchpid)
grid on 
%% I found it. man!!! small mass works for better BW --- still NOT???!!! 
m31 = 1; 
% m3 = 13.2
% m3 = 1; 
k31 = 1e6;
%k3 = 1e8; 
w31 = sqrt(k31/m31)
ccr31 = 2*sqrt(k31*m31)
zeta31 = 0.005
%zeta31 = 2.664e-3
c31 = zeta31*ccr31

Ass31 = [0, 1; -k31/m31, -c31/m31];
Bss31 = [0; 1/m31];
Css31 = [1, 0];
Dss31 = 0;
ssm31 = ss(Ass31, Bss31, Css31, Dss31)
tfm31 = tf(ssm31)
% Notch filter design for better 
tfm31 = 1e0*tfm31

figure(15)
bode(tfm31)
grid on 

figure(16)
nyquist(tfm31)
grid on 

figure(17)
rlocus(tfm31)
grid on 
% first design notch to supress the internal dynamics for stability
% purposes 
% Notch is also called: Band-Stop - Customized
wcnotch31 = w31; % Location of Notch. Notch cutoff / corner freq.
zetanotch31 = 2e-3; % Depth of Notch - bigger value ends up in blunt/dull deptth
zetanotch31 = 5e-3; 
anotch31 = 0.05; % Width of Notch - bigger value ends up in narrower filter  
anotch31 = 1e-2;
notchtfp311 = tf([1, 2*zetanotch31*wcnotch31, wcnotch31^2], [0, 0, wcnotch31^2]);
notchtfp321 = tf([0, 0, anotch31*wcnotch31], [1, anotch31*wcnotch31]);
notchtfp331 = tf([0, 0, wcnotch31/anotch31], [1, wcnotch31/anotch31]);

notchtfall31 = notchtfp311 *  notchtfp321 * notchtfp331
tfm31notch = tfm31 *  notchtfall31

figure(18)
bode(tfm31)
hold on 
bode(tfm31notch)
hold off
% %legend({'Model 3 without notch', 'Model 3 with Notch Filter - Loop-Shaping'}, 'FontSize', 14)
grid on

% tfpidm31 = pid(1e1, 5e6, 0)
tfpidm31 = pid(2045183, 0, 0)
tfm31pid = tfm31 * tfpidm31
tfm31notchpid = tfm31 *  notchtfall31 * tfpidm31
% Very important: too much strong and deep / steep notch will reduce BW
% which is BAD for ref. tracking. for stability analysis only that is okay
% but ...
% Be careful as pidTuner app shows open loop step fcn but we scrutinize the
% closed lop tf 
% Note: step response with notch is better in terms of stability and less
% oscillations but of course not good ref. tracking. for ref. tracking you
% have to do: ctrl with big kp - in pidTuner app; interestingly easily by
% increasin kp it will become unstable - now i have BW issue, bw is low !!!
% increasing Kp is also not very effective and BW does not have an explicit
% formula to do compensation design ...
% for bode oltf for step cltf
cltfm31 = feedback(tfm31, 1, -1)
cltfm31notch = feedback(tfm31notch, 1, -1)
cltfm31pid = feedback(tfm31pid, 1, -1)
cltfm31notchpid = feedback(tfm31notchpid, 1, -1)
figure(20)
step(cltfm31)
hold on 
step(cltfm31notch)
hold on 
step(cltfm31pid)
hold on
step(cltfm31notchpid)
hold off
grid on 
%legend({'Model 3 without notch', 'Model 3 with Notch Filter - Loop-Shaping', ...
    %'Model 3 without Notch but With PID', 'Model 3 with Notch and PID'}, 'FontSize', 14)
xlim([0, 10])
ylim([0, 1.25])

figure(21)
bode(cltfm31)
hold on 
bode(cltfm31notch)
hold on 
bode(cltfm31pid)
hold on
bode(cltfm31notchpid)
hold off
grid on 
%legend({'Model 3 without notch', 'Model 3 with Notch Filter - Loop-Shaping', ...
 %   'Model 3 without Notch but With PID', 'Model 3 with Notch and PID'}, 'FontSize', 14)

%% Section / project 4 - reference tracking with controller with compensator instead of notch filter for loop-shaping 
%% Section / project 3 (with ref. tracking and only collocated transfers) 
% BW is terrible idk, so go for Control1.m and look at the new design
% section and that is the project in there with good BW. 
m41 = 32;
m42 = 5;
k41 = 1e6;
k42 = 1e4; 
wn41_calc = sqrt(k41/m41)
wn42_calc = sqrt(k42/m42)
zeta41 = 5e-3;% 5e-5;
zeta42 = 25e-4;%25e-6;
ccr41 = 2*sqrt(k41*m41);
ccr42 = 2*sqrt(k42*m42);
c41 = zeta41*ccr41;
c42 = zeta42*ccr42;
m43 = 2.5; % mass of the spring added to the k2 
k43 = 1e9; % this way we shift the next res. freq. to very high values making 
% the sys as a ONE Rigid Body Mass 
mtot4 = m41+m42+m43;
A_ss42 = [0, 1;
-(k41+k42+k42)/mtot4, -(c41+c42)/mtot4] 
B_ss42 = [0;
1/mtot4];
C_ss42 = [1, 0];
C_ssnc4 = [1, 0]; % Non-collocated measurement 
D_ss42 = 0; 
wca43 = sqrt(k41/mtot4) 

ssm42 = ss(A_ss42, B_ss42, C_ss42, D_ss42) % mss stnads for state-space model

figure(41)
step(ssm42)
title('Note: this is Open loop step wo any ctrl ...')
xlim([0, 10]);
grid on 
StepInfossm42 = stepinfo(ssm42)
% No need to use pidTuner . I will tune the pid with editor command 
tfm42 = tf(ssm42)
% Design a notch for suppressing the peak at res freq. otherwise it can
% make sys unstable 
wcnotchss42 = 161; % Location of Notch. Notch cutoff / corner freq.
zetanotchss42 = 0.005; % Depth of Notch - bigger value ends up in blunt deptth
anotchss42 = 0.01; % Width of Notch - bigger value ends up in narrower filter  
notchtfp1ss42 = tf([1, 2*zetanotchss42*wcnotchss42, wcnotchss42^2], [0, 0, wcnotchss42^2]);
notchtfp2ss42 = tf([0, 0, anotchss42*wcnotchss42], [1, anotchss42*wcnotchss42]);
notchtfp3ss42 = tf([0, 0, wcnotchss42/anotchss42], [1, wcnotchss42/anotchss42]);
notchtfss42 = notchtfp1ss42*notchtfp2ss42*notchtfp3ss42;% series(notchtfp1, notchtfp2, notchtfp3); 

tfmnotchss42 = tfm42* notchtfss42;

%tfcpid42 = pid(1e4, 5e6, 1e2);
%tfcpid42 = pid(1e6, 5e6, 1e4);
%tfcpid42 = pid(1e7, 5e6, 0);
tfcpid42 = pid(1e3, 0, 0); % more than 1e7 due to sat. limits on the mot is not
% possible. and also i increased BW as much as possible by increasing the
% Kp to have good BW then design compensators to address the bad ts
oltfm42notchpid = tfm42* notchtfss42*tfcpid42; 
oltfm42pid = tfm42 * tfcpid42; 

figure(42) 
bode(ssm42) 
grid on 
hold on 
bode(tfmnotchss42) 
hold on 
bode(oltfm42notchpid) 
hold on 
bode(oltfm42pid) 
hold off 
%legend({'New design MSD', 'MSD with Notch Filter', 'MSD with Notch and P CTRL' ...
   % 'MSD with PID CTRL only'}, 'FontSize', 16)

[Gmarss42, Pmarss42, WcroGss42, WcroPss42] = margin(oltfm42pid)
% recap: do all stability analysis like on bode / nyquist / Root Locus on OL and time
% response with Mp,ess, ts ... on CL 
cltfm42notchpid = feedback(oltfm42notchpid, 1, -1);
cltfm42notch = feedback(tfmnotchss42, 1, -1);
cltfm42 = feedback(tfm42, 1, -1);
cltfm42pid = feedback(oltfm42pid, 1, -1);
figure(44)
step(cltfm42notchpid)
xlim([0, 10])
grid on 
hold on 
step(cltfm42notch)
hold on 
step(cltfm42)
hold on 
step(cltfm42pid)
%legend({'new design with Notch and P', 'new design with notch', 'new design', 'new design with p ctrl only'}, 'FontSize', 14)

SIAll = stepinfo(cltfm42pid)   % step info for all second design 
% Kp alone is terrible in ref track with significant e_ss so add I t oget
% rid of the offset 
% you can sue the pidTuner app but be very careful and import only and only
% CLTF  tfm2* notchtfss2;
% save('tfm2.mat', 'tfm2')
% save('notchtfss2.mat', 'notchtfss2')

% controlSystemDesigner(oltfm2notchpid) 
% pidTuner(oltfm2notchpid) % useless!!! unless you check for Bode only BUT NOT the stpe time 
% Does not need any compensators ... well design a compensator when you do
% not have any reference tracking ... 
% but indeed worked with pid as well so we can now design a compensator or
% better do fine-tuning on the pid 
Senss42 = 1 / (1+cltfm42pid)
figure(45)
bode(Senss42)
title("Sensitivity of 2nd Design - Accepetable Value of Sen is: mMAX of 5 dB ")
grid on 
figure(46)
rlocus(cltfm42pid)
grid on 
figure(47)
nyquist(cltfm42pid)
grid on 
% at the end decided to go for a desgin from the Ogata book 

%% project 4 / section 4 
% model of DC motor
m5 = 1; 
k5 = 10;
k5 = 0;
c5 = 1; 
ccr5 = 2*sqrt(k5*m5)
zeta5 = c5/ccr5
zeta5 = 0.1581; % by root locus or complex plane method as this is not a mass-spring-damper 
wn5 = sqrt(k5/m5)
wn5 = 3.1623; 
wd5 = wn5*sqrt(1-zeta5^2)
kp5=4;
kp5=10; 
Ass5 = [0, 1; -k5/m5, -c5/m5];
Bss5 = [0; 1/m5];
Css5 = [1, 0];
Dss5=0;
ssm5 = ss(Ass5, Bss5, Css5, Dss5)
tfm5 = tf(ssm5)
figure(511)
step(tfm5)
grid on 
att5 = zeta5*wn5; 
ts5 = 4/att5 % setling time with 2% criteria (3% criteria is: 3/att) - performanc ecriteria 
Mp5 = 100*exp((-zeta5*pi)/(sqrt(1-zeta5^2))) % overshoot criteria for stability and robustness as the complementaroy level for stability  
Dests5 = 2.6667; 
DesMp5 = 16.3034;
syms zetaG
MpG = 100*exp((-zetaG*pi)/(sqrt(1-zetaG^2))) % General formula for overshoot 
zetaF  = solve( MpG == 16.3034 , zetaG)
% if zetaF > 0 
%     zetaF = zetaF
% end
Deszeta5 = double(abs(zetaF(1,1)))
Deswn5 = 4/(Dests5*Deszeta5)
Desatt5 = Deszeta5*Deswn5
Deswd5 = Deswn5*sqrt(1-Deszeta5^2)
% from desired performance and staility criteria (desired ts and Mp; find
% the desired damping ratio and the undamped natural freq, from these find
% the desired cl poles locations, sketch on root-locus, find the angle
% deficiency which needs to be contributed by the lead C., then both params
% are determined here, how? 2 methods: geometric approach bi-sections, and
% locate zero of lead c somehow to cancel out the pole of the plant tf. afterwards check the ma
% magnitude condition to find the Kc gain of the compenasator. find the
% param of the lag c by static velocity error constnat., also do sanity check on angle and
% magnitude of the lag to be small than 5 degrees and mag around 1. and
% lastly choose the last par of the lag somehow to meet the angle condition
% and magnitude condition of the lag c. 
DesS51 = -Desatt5 + 1j*Deswd5;
DesS52 = -Desatt5 - 1j*Deswd5;
beta5 = rad2deg(acos(Deszeta5))
% angle Deficiency - thisis is exape at page 336 pgata note that the plant
%  has Kp or proportional controller as well 
tfpid5 = pid(kp5, 0, 0)
tfm5pid5 = tfm5 * tfpid5
angletfm5atDesS51 = rad2deg(angle(evalfr(tfm5pid5, DesS51))) 
angleDef5 = 180-angletfm5atDesS51
% note if angle def is > 90, you need multiple lead C 
LeCZloc5 = Deswn5 * (sin(deg2rad(beta5-angleDef5/2)) / sin(deg2rad(180-beta5-(beta5-angleDef5/2))))
LeCPloc5 = Deswn5 * (sin(deg2rad(beta5+angleDef5/2)) / sin(deg2rad(180-beta5-(beta5+angleDef5/2))))
% do magnitude condition to find the KComp
tfLeadComp = tf([1, LeCZloc5], [1, LeCPloc5])
tfolm5pidleadc = tfm5 * tfpid5 * tfLeadComp
KComp = 1 / abs((evalfr(tfolm5pidleadc, DesS51)))
% Now it's about the time to find the lag c
% static velocity error constant 
% but I go with static position error constant 
currentSPEC = (k5-kp5) / k5 % is zero 
% tfLeadComp = KComp*tf([1, LeCZloc5], [1, LeCPloc5])
tfolm5pidleadc = tfm5 * tfpid5 * tfLeadComp
tfclm5pidleadc = feedback(tfolm5pidleadc, 1, -1)
figure(512)
step(tfclm5pidleadc)
grid on 
% there is offset - steady-state error so we need lag c to address this and
% improve steady-state behavior 
tfLagComp = tf([1, 0.2], [1, 0.2*1.5])
tfolm5pidleadlagc = tfm5 * tfpid5 * tfLeadComp * tfLagComp
tfclm5pidleadlagc = feedback(tfolm5pidleadlagc, 1, -1)
figure(513)
step(tfclm5pidleadlagc)
tfclm5 = feedback(tfm5, 1, -1);
tfclm5pid5 = feedback(tfm5pid5, 1, -1)
figure(514)
step(tfclm5)
hold on 
step(tfclm5pid5)
hold on 
step(tfclm5pidleadc)
hold on
step(tfclm5pidleadlagc)
hold off
grid on
%legend({'mass-damper/DC motor', 'plant with proportional', 'plant with Kp with Lead C',...
   % 'plant with Kp with Lead-Lag C'}, 'FontSize', 14)

%%%% HOORAAAAAA it worked 
% so for mass-damper system which is equaivalent to DC motor, compensator
% with Kp proportional ctrl worked for ref tracking along with better
% performance criteria. (only difference is to find the undamped frequency and dampin
% g ratio from complex plane). 
%note lead zero and pole should be as far as possible from each other to
%change the tf as much as possible to imrpve transients and lag zero and
%pole should be very close to each other to make it less changed in
%transients but to change the offset amuch as pssible 
% only focus on the examples provided by Ogata only 
figure(515)
bode(tfm5)
hold on 
bode(tfm5pid5)
hold on 
bode(tfolm5pidleadc)
hold on
bode(tfolm5pidleadlagc)
hold off
grid on
%legend({'mass-damper/DC motor', 'plant with proportional', 'plant with Kp with Lead C',...
  %  'plant with Kp with Lead-Lag C'}, 'FontSize', 14)
% So I was able to apply lead lag c on proportional ctrl for mass-damper or
% dc motor ref tracking!!! HPOOOOOOORAAAAAAAA
% ***** Unit-Step Response of Compensated and Uncompensated Systems *****
num1 = [12.287 23.876];
den1 = [1 5.646 16.933 23.876];
num2 = [9];
den2 = [1 3 9];
num = [10];
den = [1 1 10];
t = 0:0.05:5;
c1 = step(num1,den1,t);
c2 = step(num2,den2,t);
c = step(num,den,t);
figure(55)
plot(t,c1,'-',t,c2,'.',t,c,'x')
grid on 
%legend({'Compensated System - Method: Geometry/Trigonometry', 'Compensated System - ...' ...
 %   'Method: Zero/Pole Cancelation by Lead C', 'Uncompensated System'}, 'FontSize', 12)
title('Unit-Step Responses of Compensated Systems and Uncompensated System')
xlabel('t Sec')
ylabel('Outputs c1, c2, and c')
text(1.51,1.48,'Compensated System (Method 1)')
text(0.9,0.48,'Compensated System (Method 2)')
text(2.51,0.67,'Uncompensated System')
% in sohrt; wo/ kp no ref tracking, wo/ ki no 0 offset. with both of them
% very complicated and not good idea to add comp. so go with dc motor
% models like mass-damper and no sppring, find wn and zeta rom the
% root-locus 
