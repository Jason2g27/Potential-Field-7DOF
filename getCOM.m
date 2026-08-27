function pos_cm = getCOM(mass, arm_length, cm_length, state)

    t1 =  state(1);
    t2 =  state(2);
    t3 =  state(3);
    t4 =  state(4);
    lu = arm_length(1);
    lf = arm_length(2);
    lw = arm_length(3);
    
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

    P3c3 = [luc;0;0];
    P5c5 = [lfc;0;0];
    P7c7 = [lwc;0;0];
    T0c3 = T03*[eye(3) P3c3;0 0 0 1];
    P0c3 = T0c3(1:3,4);
    T0c5 = T0p*[eye(3) P5c5;0 0 0 1];
    P0c5 = T0c5(1:3,4);
    T0c7 = T07*[eye(3) P7c7;0 0 0 1];
    P0c7 = T0c7(1:3,4);


    T3cm1 = DHConvention(0, 0, cm_length(1), 0);
    T0cm1 = T01 * T12 * T23 * T3cm1;
    T0cm2 = T01 * T12 * T23 * T34 * DHConvention(0, cm_length(2), 0, 0);
    pos_cm = (mass(1)*T0cm1(1:3,4) + mass(2)*T0cm2(1:3,4))/(sum(mass));
    
end