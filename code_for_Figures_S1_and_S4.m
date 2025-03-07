% This script takes a point in the line attractor (variable n_target), a
% range of velocities, an eta, and returns p_n and t_n. It also takes the length of the track, 
% and duration of a lap to show results in terms of range distance and
% range duration
% Employed to construct Figures S1 and S4

%% range duration and range duration for a range of velocities, for a given 
% neuron in the attractor, with a given eta. Employed to construct Figure S1
clearvars
close all
v_min=5;    %lowest velocity in the range (cm/seg)
v_max=28;   %highest velocity in the range (cm/seg)
v_pasos=100;%steps from lowest to highest velocity in the range
n_target=50;%variable n (the chosen neuron, a point in the line attractor)
eta=0.5;
l_total=150; %the length of the track (cm)
t_total=10;  %the duration of the lap (seg)

ancho=0.15;
alto=0.2;

[t_n, p_n, z_t_disparo, z_p_disparo, v, v_n]=compute_field_centres(n_target, eta, v_min, v_max, v_pasos);
r_duration = l_total./v;
r_distance = t_total*v;

fig1 = figure;
sp1=subplot(1,3,1);
plot(v, v_n)
xlabel('velocity (cm/s)')
ylabel("attractor's velocity (a.u.)")

sp2=subplot(1,3,2);
plot(p_n, r_duration)
xlabel('position (cm)')
ylabel('range duration (s)')
xlim([0, l_total])

sp3=subplot(1,3,3);
plot(t_n, r_distance)
xlabel('time (s)')
ylabel('range distance (cm)')
xlim([0, t_total])

sp1.Position=[sp1.Position(1),sp1.Position(2), ancho, alto];
sp2.Position=[sp2.Position(1),sp2.Position(2), ancho, alto];
sp3.Position=[sp3.Position(1),sp3.Position(2), ancho, alto];


%% range duration and range distance of the case of a pure place cell. Employed 
% to construct Figure S4

p_disparo_medio = mean(p_n);
t_disparo_indep = p_disparo_medio./v;
r_duration = l_total./v;
r_distance = t_total*v;

figure;

sp2=subplot(1,2,1);
plot(ones(size(r_duration))*p_disparo_medio, r_duration)
xlabel('position (cm)')
ylabel('range duration (s)')
title('posicion de disparo')    
xlim([0, l_total])

sp3=subplot(1,2,2);
plot(t_disparo_indep, r_distance)
xlabel('time (s)')
ylabel('range distance (cm)')
title('tiempo de disparo')
xlim([0, t_total])

sp1.Position=[sp1.Position(1),sp1.Position(2), ancho, alto];
sp2.Position=[sp2.Position(1),sp2.Position(2), ancho, alto];
sp3.Position=[sp3.Position(1),sp3.Position(2), ancho, alto];
%% range duration and range distance of the case of a pure time cell. Employed 
% to construct Figure S4

t_disparo_medio = mean(t_n);
p_disparo_t_correct = v*t_disparo_medio;
r_duration = l_total./v;
r_distance = t_total*v;

fig4 = figure;


sp2=subplot(1,2,1);
plot(p_disparo_t_correct, r_duration)
xlabel('space (cm)')
ylabel('range duration (seg)')
title('posicion de disparo')    
xlim([0, l_total])

sp3=subplot(1,2,2);
plot(ones(size(r_distance))*t_disparo_medio, r_distance)
xlabel('time (seg)')
ylabel('range distance (cm)')
title('tiempo de disparo')
xlim([0, t_total])

sp1.Position=[sp1.Position(1),sp1.Position(2), ancho, alto];
sp2.Position=[sp2.Position(1),sp2.Position(2), ancho, alto];
sp3.Position=[sp3.Position(1),sp3.Position(2), ancho, alto];

    
