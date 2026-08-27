
syms t1 t2 t3 t4 t5 t6 t7 lu lf luc lfc lwc PI real

T0_p1 =[ cos(t1), cos(t2 + t3)*sin(t1), -sin(t2 + t3)*sin(t1), 0
       0,        -sin(t2 + t3),         -cos(t2 + t3), 0
-sin(t1), cos(t2 + t3)*cos(t1), -sin(t2 + t3)*cos(t1), 0
       0,                    0,                     0, 1];
 
 
T0_p2 =
 
[cos(t2 + t3)*sin(t1), sin(t2 + t3)*sin(t4 + t5)*sin(t1) - cos(t4 + t5)*cos(t1), - sin(t4 + t5)*cos(t1) - cos(t4 + t5)*sin(t2 + t3)*sin(t1),  lu*cos(t1)
       -sin(t2 + t3),                                cos(t2 + t3)*sin(t4 + t5),                                 -cos(t2 + t3)*cos(t4 + t5),           0
cos(t2 + t3)*cos(t1), cos(t4 + t5)*sin(t1) + sin(t2 + t3)*sin(t4 + t5)*cos(t1),   sin(t4 + t5)*sin(t1) - cos(t4 + t5)*sin(t2 + t3)*cos(t1), -lu*sin(t1)
                   0,                                                        0,                                                          0,           1];

