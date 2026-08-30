Te1= [-cos(t2)*sin(t1), cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2),   cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3), 0
         sin(t2),                           cos(t2)*cos(t3),                            -cos(t2)*sin(t3), 0 %Equation 4.6
-cos(t1)*cos(t2), cos(t1)*cos(t3)*sin(t2) - sin(t1)*sin(t3), - cos(t3)*sin(t1) - cos(t1)*sin(t2)*sin(t3), 0
               0,                                         0,                                           0, 1]; 


Y1 = [-Lu*cos(t2)*sin(t1);
         Lu*sin(t2); %Equation 4.7
-Lu*cos(t1)*cos(t2);
                  1];  
 
 
Y2 = [- Lf*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) - Lu*cos(t2)*sin(t1);
                                              Lf*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)) + Lu*sin(t2); %Equation 4.8
  Lf*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - Lu*cos(t1)*cos(t2);
                                                                                                        1];
 
 
 
Tdy1 = [Lw*(sin(t5)*sin(t7) - cos(t5)*cos(t7)*sin(t6));
                     - Lf - Lw*cos(t6)*cos(t7); %Equation 4.11
Lw*(cos(t5)*sin(t7) + cos(t7)*sin(t5)*sin(t6));
                                             1];
 
 
Tdy2 = [Lw*(sin(t5)*sin(t7) - cos(t5)*cos(t7)*sin(t6)) + Lt*(cos(t8)*(sin(t5)*sin(t7) - cos(t5)*cos(t7)*sin(t6)) + sin(t8)*(cos(t7)*sin(t5) + cos(t5)*sin(t6)*sin(t7)));
                                                                               Lt*(cos(t6)*sin(t7)*sin(t8) - cos(t6)*cos(t7)*cos(t8)) - Lf - Lw*cos(t6)*cos(t7); %Equation 4.12
Lw*(cos(t5)*sin(t7) + cos(t7)*sin(t5)*sin(t6)) + Lt*(cos(t8)*(cos(t5)*sin(t7) + cos(t7)*sin(t5)*sin(t6)) + sin(t8)*(cos(t5)*cos(t7) - sin(t5)*sin(t6)*sin(t7)));
                                                                                                                                                              1];
 
 
Tdz = [cos(t2)*sin(t1)*sin(t4) - cos(t1)*cos(t3)*cos(t4) + cos(t4)*sin(t1)*sin(t2)*sin(t3),   cos(t2)*cos(t4)*sin(t3) - sin(t2)*sin(t4), cos(t1)*cos(t2)*sin(t4) + cos(t3)*cos(t4)*sin(t1) + cos(t1)*cos(t4)*sin(t2)*sin(t3), Lu*sin(t4);
cos(t2)*cos(t4)*sin(t1) + cos(t1)*cos(t3)*sin(t4) - sin(t1)*sin(t2)*sin(t3)*sin(t4), - cos(t4)*sin(t2) - cos(t2)*sin(t3)*sin(t4), cos(t1)*cos(t2)*cos(t4) - cos(t3)*sin(t1)*sin(t4) - cos(t1)*sin(t2)*sin(t3)*sin(t4), Lu*cos(t4);   %Equation 4.14   
                                          cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2),                             cos(t2)*cos(t3),                                           cos(t1)*cos(t3)*sin(t2) - sin(t1)*sin(t3),          0;
                                                                                  0,                                           0,                                                                                   0,          1];
 
                                                                                                                                                         