function G = getG(model, state)
    lu = model.lu;
    lf = model.lf;
    lw = model.lw;
    luc = model.Lcmu;
    lfc = model.Lcmf;
    lwc = model.Lcmw;
    t1 = state(1);
    t2 = state(2);
    t3 = state(3);
    t4 = state(4);
    t5 = state(5);
    t6 = state(6);
    t7 = state(7);
    td1 = state(8);
    td2 = state(9);
    td3 = state(10);
    td4 = state(11);
    td5 = state(12);
    td6 = state(13);
    td7 = state(14);
    Ixx_3 = model.Ic33(1,1);
    Iyy_3 = model.Ic33(2,2);
    Izz_3 = model.Ic33(3,3);
    Ixx_5 = model.Ic55(1,1);
    Iyy_5 = model.Ic55(2,2);
    Izz_5 = model.Ic55(3,3);
    Ixx_7 = model.Ic77(1,1);
    Iyy_7 = model.Ic77(2,2);
    Izz_7 = model.Ic77(3,3);
    m3 = model.m3;
    m5 = model.m5;
    m7 = model.m7;
    g = model.g;


G=[g*lu*m5*cos(t2)*sin(t1) + g*lu*m7*cos(t2)*sin(t1) + g*luc*m3*cos(t2)*sin(t1) + g*lf*m7*cos(t2)*cos(t4)*sin(t1) + g*lf*m7*cos(t1)*cos(t3)*sin(t4) + g*lfc*m5*cos(t2)*cos(t4)*sin(t1) + g*lfc*m5*cos(t1)*cos(t3)*sin(t4) - g*lwc*m7*cos(t1)*cos(t5)*sin(t3)*sin(t7) - g*lf*m7*sin(t1)*sin(t2)*sin(t3)*sin(t4) - g*lfc*m5*sin(t1)*sin(t2)*sin(t3)*sin(t4) - g*lwc*m7*cos(t3)*cos(t5)*sin(t1)*sin(t2)*sin(t7) - g*lwc*m7*cos(t1)*cos(t7)*sin(t3)*sin(t5)*sin(t6) - g*lwc*m7*cos(t2)*sin(t1)*sin(t4)*sin(t5)*sin(t7) + g*lwc*m7*cos(t2)*cos(t4)*cos(t6)*cos(t7)*sin(t1) + g*lwc*m7*cos(t1)*cos(t3)*cos(t6)*cos(t7)*sin(t4) + g*lwc*m7*cos(t1)*cos(t3)*cos(t4)*sin(t5)*sin(t7) - g*lwc*m7*cos(t1)*cos(t3)*cos(t4)*cos(t5)*cos(t7)*sin(t6) + g*lwc*m7*cos(t2)*cos(t5)*cos(t7)*sin(t1)*sin(t4)*sin(t6) - g*lwc*m7*cos(t6)*cos(t7)*sin(t1)*sin(t2)*sin(t3)*sin(t4) - g*lwc*m7*cos(t3)*cos(t7)*sin(t1)*sin(t2)*sin(t5)*sin(t6) - g*lwc*m7*cos(t4)*sin(t1)*sin(t2)*sin(t3)*sin(t5)*sin(t7) + g*lwc*m7*cos(t4)*cos(t5)*cos(t7)*sin(t1)*sin(t2)*sin(t3)*sin(t6)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              g*cos(t1)*(lu*m5*sin(t2) + lu*m7*sin(t2) + luc*m3*sin(t2) + lf*m7*cos(t4)*sin(t2) + lfc*m5*cos(t4)*sin(t2) + lf*m7*cos(t2)*sin(t3)*sin(t4) + lfc*m5*cos(t2)*sin(t3)*sin(t4) + lwc*m7*cos(t2)*cos(t3)*cos(t5)*sin(t7) + lwc*m7*cos(t4)*cos(t6)*cos(t7)*sin(t2) - lwc*m7*sin(t2)*sin(t4)*sin(t5)*sin(t7) + lwc*m7*cos(t2)*cos(t6)*cos(t7)*sin(t3)*sin(t4) + lwc*m7*cos(t2)*cos(t3)*cos(t7)*sin(t5)*sin(t6) + lwc*m7*cos(t2)*cos(t4)*sin(t3)*sin(t5)*sin(t7) + lwc*m7*cos(t5)*cos(t7)*sin(t2)*sin(t4)*sin(t6) - lwc*m7*cos(t2)*cos(t4)*cos(t5)*cos(t7)*sin(t3)*sin(t6))
                                                                                                                                                                                                                                                                                                                     g*lf*m7*cos(t1)*cos(t3)*sin(t2)*sin(t4) - g*lfc*m5*sin(t1)*sin(t3)*sin(t4) - g*lf*m7*sin(t1)*sin(t3)*sin(t4) + g*lfc*m5*cos(t1)*cos(t3)*sin(t2)*sin(t4) - g*lwc*m7*cos(t3)*cos(t5)*sin(t1)*sin(t7) - g*lwc*m7*cos(t1)*cos(t5)*sin(t2)*sin(t3)*sin(t7) - g*lwc*m7*cos(t6)*cos(t7)*sin(t1)*sin(t3)*sin(t4) - g*lwc*m7*cos(t3)*cos(t7)*sin(t1)*sin(t5)*sin(t6) - g*lwc*m7*cos(t4)*sin(t1)*sin(t3)*sin(t5)*sin(t7) + g*lwc*m7*cos(t1)*cos(t3)*cos(t6)*cos(t7)*sin(t2)*sin(t4) + g*lwc*m7*cos(t1)*cos(t3)*cos(t4)*sin(t2)*sin(t5)*sin(t7) + g*lwc*m7*cos(t4)*cos(t5)*cos(t7)*sin(t1)*sin(t3)*sin(t6) - g*lwc*m7*cos(t1)*cos(t7)*sin(t2)*sin(t3)*sin(t5)*sin(t6) - g*lwc*m7*cos(t1)*cos(t3)*cos(t4)*cos(t5)*cos(t7)*sin(t2)*sin(t6)
                                                                                                                                                                                                                                                                                           g*lf*m7*cos(t1)*cos(t2)*sin(t4) + g*lf*m7*cos(t3)*cos(t4)*sin(t1) + g*lfc*m5*cos(t1)*cos(t2)*sin(t4) + g*lfc*m5*cos(t3)*cos(t4)*sin(t1) + g*lf*m7*cos(t1)*cos(t4)*sin(t2)*sin(t3) + g*lfc*m5*cos(t1)*cos(t4)*sin(t2)*sin(t3) - g*lwc*m7*cos(t3)*sin(t1)*sin(t4)*sin(t5)*sin(t7) + g*lwc*m7*cos(t1)*cos(t2)*cos(t6)*cos(t7)*sin(t4) + g*lwc*m7*cos(t3)*cos(t4)*cos(t6)*cos(t7)*sin(t1) + g*lwc*m7*cos(t1)*cos(t2)*cos(t4)*sin(t5)*sin(t7) - g*lwc*m7*cos(t1)*cos(t2)*cos(t4)*cos(t5)*cos(t7)*sin(t6) + g*lwc*m7*cos(t1)*cos(t4)*cos(t6)*cos(t7)*sin(t2)*sin(t3) + g*lwc*m7*cos(t3)*cos(t5)*cos(t7)*sin(t1)*sin(t4)*sin(t6) - g*lwc*m7*cos(t1)*sin(t2)*sin(t3)*sin(t4)*sin(t5)*sin(t7) + g*lwc*m7*cos(t1)*cos(t5)*cos(t7)*sin(t2)*sin(t3)*sin(t4)*sin(t6)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              g*lwc*m7*(sin(t1)*sin(t3)*sin(t5)*sin(t7) + cos(t1)*cos(t2)*cos(t5)*sin(t4)*sin(t7) + cos(t3)*cos(t4)*cos(t5)*sin(t1)*sin(t7) - cos(t1)*cos(t3)*sin(t2)*sin(t5)*sin(t7) - cos(t5)*cos(t7)*sin(t1)*sin(t3)*sin(t6) + cos(t1)*cos(t3)*cos(t5)*cos(t7)*sin(t2)*sin(t6) + cos(t1)*cos(t4)*cos(t5)*sin(t2)*sin(t3)*sin(t7) + cos(t1)*cos(t2)*cos(t7)*sin(t4)*sin(t5)*sin(t6) + cos(t3)*cos(t4)*cos(t7)*sin(t1)*sin(t5)*sin(t6) + cos(t1)*cos(t4)*cos(t7)*sin(t2)*sin(t3)*sin(t5)*sin(t6))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 -g*lwc*m7*cos(t7)*(cos(t3)*sin(t1)*sin(t4)*sin(t6) - cos(t1)*cos(t2)*cos(t4)*sin(t6) + cos(t6)*sin(t1)*sin(t3)*sin(t5) + cos(t1)*cos(t2)*cos(t5)*cos(t6)*sin(t4) + cos(t3)*cos(t4)*cos(t5)*cos(t6)*sin(t1) - cos(t1)*cos(t3)*cos(t6)*sin(t2)*sin(t5) + cos(t1)*sin(t2)*sin(t3)*sin(t4)*sin(t6) + cos(t1)*cos(t4)*cos(t5)*cos(t6)*sin(t2)*sin(t3))
                                                                                                                                                                                                                                                                                                                                                                                                                        g*lwc*m7*(cos(t1)*cos(t3)*cos(t5)*cos(t7)*sin(t2) - cos(t5)*cos(t7)*sin(t1)*sin(t3) + cos(t1)*cos(t2)*cos(t4)*cos(t6)*sin(t7) + cos(t1)*cos(t2)*cos(t7)*sin(t4)*sin(t5) + cos(t3)*cos(t4)*cos(t7)*sin(t1)*sin(t5) - cos(t3)*cos(t6)*sin(t1)*sin(t4)*sin(t7) + sin(t1)*sin(t3)*sin(t5)*sin(t6)*sin(t7) + cos(t1)*cos(t4)*cos(t7)*sin(t2)*sin(t3)*sin(t5) + cos(t1)*cos(t2)*cos(t5)*sin(t4)*sin(t6)*sin(t7) + cos(t3)*cos(t4)*cos(t5)*sin(t1)*sin(t6)*sin(t7) - cos(t1)*cos(t6)*sin(t2)*sin(t3)*sin(t4)*sin(t7) - cos(t1)*cos(t3)*sin(t2)*sin(t5)*sin(t6)*sin(t7) + cos(t1)*cos(t4)*cos(t5)*sin(t2)*sin(t3)*sin(t6)*sin(t7))];
 
end