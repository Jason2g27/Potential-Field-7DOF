 clear
clc
close all

%% For the complex case

syms q1 q2 q3 q4 q5 q6 q7 lu lf lw PI q1d q2d q3d q4d Ft Fn real



T01 = DHConvention(-PI/2, 0, 0, -PI/2 + q1); 
T01 = simplify(subs(T01, [cos(PI/2), sin(PI/2)], [0, 1]));

T12 = DHConvention(-PI/2, 0, 0, -PI/2 + q2); 
T12 = simplify(subs(T12, [cos(PI/2), sin(PI/2)], [0, 1]));

T23 = DHConvention(PI/2, 0, 0, q3); %
T23 = simplify(subs(T23, [cos(PI/2), sin(PI/2)], [0, 1]));

T3_p1= DHConvention(-PI/2, 0, 0, -PI/2);
T3_p1 = simplify(subs(T3_p1, [cos(PI/2), sin(PI/2)], [0, 1]));

Tp1_4 = DHConvention(-PI/2, lu, 0, PI/2+q4); 
Tp1_4 = simplify(subs(Tp1_4, [cos(PI/2), sin(PI/2)], [0, 1]));

T45 = DHConvention(-PI/2, 0, 0, q5); 
T45 = simplify(subs(T45, [cos(PI/2), sin(PI/2)], [0, 1]));

T5_p2 = DHConvention(-PI/2, 0, 0, PI/2); 
T5_p2 = simplify(subs(T5_p2, [cos(PI/2), sin(PI/2)], [0, 1]));

Tp2_6 = DHConvention(0, lf, 0, q6); 
Tp2_6 = simplify(subs(Tp2_6, [cos(PI/2), sin(PI/2)], [0, 1]));


T67 = DHConvention(-PI/2, 0, 0, q7); 
T67 = simplify(subs(T67, [cos(PI/2), sin(PI/2)], [0, 1]));



T78= DHConvention(0, lw, 0, 0); 
T78 = simplify(subs(T78, [cos(PI/2), sin(PI/2)], [0, 1]));


T02 = T01*T12;
T02 = simplify(subs(T02, [cos(PI/2), sin(PI/2)], [0, 1]));
T02 = simplify(subs(T02, PI, PI));

T03 = T02*T23;
T03 = simplify(subs(T03, [cos(PI/2), sin(PI/2)], [0, 1]));
T03 = simplify(subs(T03, PI, PI));

T0_p1= T03 * T3_p1;
T04 = T0_p1 * Tp1_4;
T04 = simplify(subs(T04, [cos(PI/2), sin(PI/2)], [0, 1]));
T04 = simplify(subs(T04, PI, PI));

T05 = T04 * T45;
T05 = simplify(subs(T05, [cos(PI/2), sin(PI/2)], [0, 1]));
T05 = simplify(subs(T05, PI, PI));

T0_p2 = T05 * T5_p2;
T06 = T0_p2 * Tp2_6;
T06 = simplify(subs(T06, [cos(PI/2), sin(PI/2)], [0, 1]));
T06 = simplify(subs(T06, PI, PI));

T07 = T06 * T67;
T07 = simplify(subs(T07, [cos(PI/2), sin(PI/2)], [0, 1]));
T07 = simplify(subs(T07, PI, PI));

T08=T07*T78;
T08 = simplify(subs(T08, [cos(PI/2), sin(PI/2)], [0, 1]));
T08 = simplify(subs(T08, PI, PI));




Jve = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,T08,8);    
Jve = simplify(subs(Jve, PI, PI)); 

% % rank_Je = rank(simplify(Jve'*Jve));
% % for elbow contact
% Jvel = calculateJacobian(T01,T02,T03,T04,T05,4);
% Jvel = simplify(subs(Jvel, PI, PI));
% 
% v_elbow = Jvel * [q1d;q2d;q3d];
% % syms Fev Feh;
% % tau_elbow_v = simplify(Jvel'*Fev*[0;0;1])
% % 
% % tau_elbow_h = simplify(Jvel'*Feh*cross(v_elbow, [0;0;1])/norm(cross(v_elbow, [0;0;1])))
% 
syms f1 f2 f3;
% tau_elbow = simplify(Jvel'*[f1;f2;f3])
% 
% % for tangential motion
% ve = Jve *[q1d;q2d;q3d;q4d];
% z4 = T04(1:3,3)/norm(T04(1:3,3));
% % z4 and ve will never align with each other since the elbow will always
% % move and there will be motion component perpendicular to z4.
% nt = cross(ve,z4)/norm(cross(ve,z4));
% Fnt = Ft*nt;  % force in the plane
% tangentialTorque = simplify(Jve'*Fnt)
% nn = z4;
% Fnn = Fn*nn;  % force normal to the plane
% normalTorque = simplify(Jve'*Fnn)
tau_wrist = simplify(Jve'*[f1;f2;f3])
% for sanity test
% tangentialTorqueSanityTest = simplify(subs(tangentialTorque, [q1, q3, q1d, q3d], [0,0,0,0]))
% normalTorqueSanityTest = simplify(subs(normalTorque, [q1, q3, q1d, q3d], [0,0,0,0]))



%% For the simple case
% syms F q1 q2 q1d q2d lu lf real
% J = [0 0
%      lu*cos(q1)+lf*cos(q1+q2) lf*cos(q1+q2)
%      lu*sin(q1)+lf*sin(q1+q2) lf*sin(q1+q2)];
% 
% v = J*[q1d;q2d]/norm(J*[q1d;q2d]);
% z2 = [1;0;0];
% n = cross(z2,v);
% Fn = F*n;
% torq = simplify(J'*Fn)





%% Plot
% % complex version
% lf = 1;
% lu = 1.2;
% q1 = 0;
% q2 = 1.2;
% q3 = 0;
% q4 = 0.5;
% q2d = 1;
% q4d = 2;
% q3d = 0;
% q1d = 0;
% F = 1;
% Ft = F;
% Fn = F;
% T01 = DHConvention(-PI/2, 0, 0, -PI/2 + q1);
% T12 = DHConvention(-PI/2, 0, 0, -PI/2 + q2);
% T23 = DHConvention(PI/2, 0, 0, q3);
% T33p = DHConvention(-PI/2, 0, 0, -PI/2);
% T3p4 = DHConvention(0, lu, 0, q4);
% T45 = DHConvention(0, lf, 0, 0);
% T02 = T01 * T12;
% T03 = T02 * T23;
% 
% T04 = T03 * T33p * T3p4;
% 
% T05 = T04 * T45;
% Pos = [zeros(3,1) T01(1:3,4) T02(1:3,4) T03(1:3,4) T04(1:3,4) T05(1:3,4)];
% 
% Jve = calculateJacobian(T01,T02,T03,T04,T05,5);    
% 
% % for tangential motion
% ve = Jve *[q1d;q2d;q3d;q4d]/norm(Jve *[q1d;q2d;q3d;q4d]);
% z4 = T04(1:3,3)/norm(T04(1:3,3));
% nt = cross(ve,z4);
% Fnt = Ft*nt;  % force in the plane
% Fn_comp = Fnt;
% tangentialTorque = simplify(Jve'*Fnt)
% nn = cross(ve,nt);
% Fnn = Fn*nn;  % force normal to the plane
% normalTorque = simplify(Jve'*Fnn)
% torq_comp = normalTorque;
% 
% figure()
% subplot(1,2,1)
% plot3(Pos(1,:),Pos(2,:),Pos(3,:),"-o")
% hold on
% P = [Pos(:,end) Pos(:,end)+Fn];
% plot3(P(1,:), P(2,:), P(3,:), 'LineWidth', 2);
% axis equal
% xlabel("x")
% ylabel("y")
% zlabel("z")
% axis([-5 5 -5 5 -5 5])
% view(120,20)
% 
% % simple version
% q1 = q2;
% q2 = q4;
% q1d = q2d;
% q2d = q4d;
% J = [0 0
%      lu*cos(q1)+lf*cos(q1+q2) lf*cos(q1+q2)
%      lu*sin(q1)+lf*sin(q1+q2) lf*sin(q1+q2)];
% 
% v = J*[q1d;q2d]/norm(J*[q1d;q2d]);
% pe = [0;lu*sin(q1)+lf*sin(q1+q2);-lu*cos(q1)-lf*cos(q1+q2)]/norm([0;lu*sin(q1)+lf*sin(q1+q2);-lu*cos(q1)-lf*cos(q1+q2)]);
% n = (pe-dot(pe,v)*v)/norm(pe-dot(pe,v)*v);
% % n = cross(v,[1;0;0]);
% Fn = F*n;
% torq = J'*Fn;
% Pos = [zeros(3,1) [0;lu*sin(q1);-lu*cos(q1)] [0;lu*sin(q1)+lf*sin(q1+q2);-lu*cos(q1)-lf*cos(q1+q2)]];
% subplot(1,2,2)
% plot3(Pos(1,:),Pos(2,:),Pos(3,:),"-o")
% hold on
% P = [Pos(:,end) Pos(:,end)+Fn];
% plot3(P(1,:), P(2,:), P(3,:), 'LineWidth', 2);
% axis equal
% xlabel("x")
% ylabel("y")
% zlabel("z")
% axis([-5 5 -5 5 -5 5])
% view(120,20)
% difference = norm(Fn_comp - Fn);
% fprintf('The difference between two reaction forces is %d\n', difference);
% difference = norm(torq_comp - [0; torq(1); 0 ;torq(2)]);
% fprintf('The difference between two torques is %d\n', difference);




