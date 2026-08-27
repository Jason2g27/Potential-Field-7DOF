% Arm dynamics

clear
clc
close all
syms t1 t2 t3 t4 t5 t6 t7 lu lf lw luc lfc lwc PI real
PI=sym(pi);
T01 = DHConvention(-PI/2, 0, 0, -PI + t1); %
T12 = DHConvention(-PI/2, 0, 0, -PI/2 + t2); %
T02 = T01*T12;
T23 = DHConvention(PI/2, 0, 0, -PI/2+t3); %
T03 = T02*T23;
T34 = DHConvention(-PI/2, lu, 0, PI/2+t4); %
T04 = T03 * T34;
T45 = DHConvention(-PI/2, 0, 0, t5); 
T05 = T04 * T45;
T5p = DHConvention(-PI/2, 0, 0, PI/2); 
T0p = T05 * T5p;
Tp6 = DHConvention(0, lf, 0, t6); 
T06 = T0p * Tp6;
T67 = DHConvention(-PI/2, 0, 0, t7); 
T07 = T06 * T67;
T78= DHConvention(0, lw, 0, 0); 
T08=T07*T78;

%% Dynamics

syms m3 m7 m5 g real
syms Ixx_3 Iyy_3 Izz_3 Ixx_7 Iyy_7 Izz_7 Ixx_5  Iyy_5 Izz_5 real
syms td1 tdd1 td2 tdd2 td3 tdd3 td4 tdd4 tdd5 td6 tdd6 td7 tdd7 real

%assume(m3, m7, m5, lu, lf, luc, lfc, lwc,'positive')

%% Inertia
I3c3 = [Ixx_3 0 0;0 Iyy_3 0;0 0 Izz_3];
I5c5 = [Ixx_5 0 0;0 Iyy_5 0;0 0 Izz_5];
I7c7 = [Ixx_7 0 0;0 Iyy_7 0;0 0 Izz_7];

I03 = T03(1:3,1:3)*I3c3*T03(1:3,1:3)';
I05 = T0p(1:3,1:3)*I5c5*T0p(1:3,1:3)';
I07 = T07(1:3,1:3)*I7c7*T07(1:3,1:3)';

P01 = T01(1:3,4);
P02 = T02(1:3,4);
P03 = T03(1:3,4);
P04 = T04(1:3,4);
P05 = T05(1:3,4);
P06 = T06(1:3,4);
P07 = T07(1:3,4);



P3c3 = [luc;0;0];
P5c5 = [lfc;0;0];
P7c7 = [lwc;0;0];
T0c3 = T03*[eye(3) P3c3;0 0 0 1];
P0c3 = T0c3(1:3,4);
T0c5 = T0p*[eye(3) P5c5;0 0 0 1];
P0c5 = T0c5(1:3,4);
T0c7 = T07*[eye(3) P7c7;0 0 0 1];
P0c7 = T0c7(1:3,4);

%%
J0v3 = [cross(T01(1:3,3),P0c3-P01) cross(T02(1:3,3),P0c3-P02) cross(T03(1:3,3),P0c3-P03) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)];
J0w3 = [T01(1:3,3) T02(1:3,3) T03(1:3,3) zeros(3,1) zeros(3,1) zeros(3,1) zeros(3,1)];

J0v5 = [cross(T01(1:3,3),P0c5-P01) cross(T02(1:3,3),P0c5-P02) cross(T03(1:3,3),P0c5-P03) cross(T04(1:3,3),P0c5-P04) cross(T05(1:3,3),P0c5-P05) zeros(3,1) zeros(3,1)];
J0w5 = [T01(1:3,3) T02(1:3,3) T03(1:3,3) T04(1:3,3) T05(1:3,3) zeros(3,1) zeros(3,1)];

J0v7 = [cross(T01(1:3,3),P0c7-P01) cross(T02(1:3,3),P0c7-P02) cross(T03(1:3,3),P0c7-P03) cross(T04(1:3,3),P0c7-P04) cross(T05(1:3,3),P0c7-P05) cross(T06(1:3,3),P0c7-P06) cross(T07(1:3,3),P0c7-P07)];
J0w7 = [T01(1:3,3) T02(1:3,3) T03(1:3,3) T04(1:3,3) T05(1:3,3) T06(1:3,3) T07(1:3,3)];

M = J0v3'*m3*J0v3+J0w3'*I03*J0w3+J0v5'*m5*J0v5+J0w5'*I05*J0w5+J0v7'*m7*J0v7+J0w7'*I07*J0w7 ;

M=simplify(M)
V1 = 0; V2 = 0; V3 = 0;V4 = 0; V5 = 0; V6 = 0; V7 = 0;
theta = [t1 t2 t3 t4 t5 t6 t7];
theta_dot = [td1 td2 td3 td4 td5 td6 td7];
for j = 1:7
    for k = 1:7
        V1 = V1 + (diff(M(1,j),theta(k))-1/2*diff(M(j,k),theta(1)))*theta_dot(k)*theta_dot(j);
        V2 = V2 + (diff(M(2,j),theta(k))-1/2*diff(M(j,k),theta(2)))*theta_dot(k)*theta_dot(j);
        V3 = V3 + (diff(M(3,j),theta(k))-1/2*diff(M(j,k),theta(3)))*theta_dot(k)*theta_dot(j);
        V4 = V4 + (diff(M(4,j),theta(k))-1/2*diff(M(j,k),theta(4)))*theta_dot(k)*theta_dot(j);
        V5 = V5 + (diff(M(5,j),theta(k))-1/2*diff(M(j,k),theta(5)))*theta_dot(k)*theta_dot(j);
        V6 = V6 + (diff(M(6,j),theta(k))-1/2*diff(M(j,k),theta(6)))*theta_dot(k)*theta_dot(j);
        V7 = V7 + (diff(M(7,j),theta(k))-1/2*diff(M(j,k),theta(7)))*theta_dot(k)*theta_dot(j);
    end
end

C = ([V1;V2;V3;V4;V5;V6;V7]);
C=simplify(C)
G1 = -m3*[0 0 -g]*J0v3(:,1) -m5*[0 0 -g]*J0v5(:,1) -m7*[0 0 -g]*J0v7(:,1);
G2 = -m3*[0 0 -g]*J0v3(:,2) -m5*[0 0 -g]*J0v5(:,2) -m7*[0 0 -g]*J0v7(:,2);
G3 = -m3*[0 0 -g]*J0v3(:,3) -m5*[0 0 -g]*J0v5(:,3) -m7*[0 0 -g]*J0v7(:,3);
G4 = -m3*[0 0 -g]*J0v3(:,4) -m5*[0 0 -g]*J0v5(:,4) -m7*[0 0 -g]*J0v7(:,4);
G5 = -m3*[0 0 -g]*J0v3(:,5) -m5*[0 0 -g]*J0v5(:,5) -m7*[0 0 -g]*J0v7(:,5);
G6 = -m3*[0 0 -g]*J0v3(:,6) -m5*[0 0 -g]*J0v5(:,6) -m7*[0 0 -g]*J0v7(:,6);
G7 = -m3*[0 0 -g]*J0v3(:,7) -m5*[0 0 -g]*J0v5(:,7) -m7*[0 0 -g]*J0v7(:,7);
G = subs([G1;G2;G3;G4;G5;G6;G7],[cos(PI/2) sin(PI/2) cos(PI/2-t1) sin(PI/2-t1) cos(PI/2-t2) sin(PI/2-t2) cos(PI/2-t3) sin(PI/2-t3) cos(PI/2-t4) sin(PI/2-t4) cos(PI/2-t5) sin(PI/2-t5) cos(PI/2-t6) sin(PI/2-t6) cos(PI/2-t7) sin(PI/2-t7)],[0 1 sin(t1) cos(t1) sin(t2) cos(t2) sin(t3) cos(t3) sin(t4) cos(t4) sin(t5) cos(t5) sin(t6) cos(t6) sin(t7) cos(t7)]);
G=simplify(G)

% tao_ext = J_geom'*[T04(1:3,1:3) eye(3);eye(3) T04(1:3,1:3)]*[0;0;0;0;0;0];
% Q1 = simplify(M(1,1)*tdd1 + M(1,2)*tdd2 + M(1,3)*tdd3 + V1 + G1) - tao_ext(1)
% Q2 = simplify(M(2,1)*tdd1 + M(2,2)*tdd2 + M(2,3)*tdd3 + V2 + G2) - tao_ext(2)
% Q3 = simplify(M(3,1)*tdd1 + M(3,2)*tdd2 + M(3,3)*tdd3 + V3 + G3) - tao_ext(3)




