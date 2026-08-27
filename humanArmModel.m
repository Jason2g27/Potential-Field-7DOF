% human arm model
clc
clear
close all
lu = 0.3;
lf = 0.33;
Lcmu = .16;
Lcmf = .17;
m3 = 0.17;
m4 = 0.18;
Ic33  = diag([0.02 0.02 0.02]);
Ic44  = diag([0.02 0.02 0.02]);

L(1) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', -pi/2, ...
        'modified', 'm', 0, 'I', zeros(3), 'r',  [0;0;0]);
L(2) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', -pi/2, ...
    'modified', 'm', 0, 'I', zeros(3), 'r', [0;0;0]);
L(3) = Revolute('d', 0, 'a', 0, 'alpha', pi/2, 'offset', 0, ...
    'modified', 'm', 0, 'I', zeros(3), 'r', [0;0;0]);
L(4) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', -pi/2, ...
    'modified', 'm', 0, 'I', zeros(3), 'r', [0;0;0]);
L(5) = Revolute('d', 0, 'a', lu, 'alpha', -pi/2, 'offset', pi/2, ...
    'modified', 'm', m3, 'I', Ic33, 'r', [Lcmu;0;0]);
L(6) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0, ...
    'modified', 'm', m4, 'I', Ic44, 'r', [Lcmf;0;0]);

L(7) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', pi/2, ...
    'modified', 'm', 0, 'I', zeros(3), 'r', [0;0;0]);
L(8) = Revolute('d', 0, 'a', lf, 'alpha', 0, 'offset', 0, ...
    'modified', 'm', m4, 'I', Ic44, 'r', [Lcmf;0;0]);

L(9) = Revolute('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0, ...
    'modified', 'm', m4, 'I', Ic44, 'r', [0;0;0]);

L(10) = Revolute('d', 0, 'a', Lcmf, 'alpha', 0, 'offset', 0, ...
    'modified', 'm', m4/3, 'I', Ic44/9, 'r', [Lcmf/5;0;0]);





threelink = SerialLink(L, 'name', '9-link Arm');

disp(threelink.n);

q = [0 0 0 0 0 0 0 0 0 0] ; 
threelink.plot(q);
threelink.teach();
%ui_elements = findobj(gcf, '-property', 'FontSize');
%set(ui_elements, 'FontName', 'Arial', 'FontSize', 10);


% threelink.gravity = [0; -9.8; 0];
% qc = zeros(6,length(state(1,:)));
% qc(1:3,:) = state(1:3,:);
% qc(5,:) = state(4,:);
% 
% qcdot = zeros(6,length(state(1,:)));
% qcdot(1:3,:) = state(5:7,:);
% qcdot(5,:) = state(8,:);
% 
% qcddot = zeros(6,length(state(1,:)));
% qcddot(1:3,:) = state(9:11,:);
% qcddot(5,:) = state(12,:);
% 
% torque_tot = threelink.rne(qc', qcdot', qcddot');
% torque_tot = torque_tot';
% ind = [1,2,3,5];
% torque = torque_tot(ind,:);