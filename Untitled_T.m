clc;
clear all;
close all;
syms t1 t2 t3 t4 t5 t6 t7 lf lu luc lw lfc lwc PI Lu Lf real
PI=sym(pi); 

        T01p = DHConvention(-PI/2, 0, 0, -PI/2 + t1); %
        T12p = DHConvention(-PI/2, 0, 0, -PI/2 + t2); %
        T02p = T01p*T12p;
        T23p = DHConvention(PI/2, 0, 0, t3); %
        T03p = T02p*T23p;
        T3_p1p= DHConvention(-PI/2, 0, 0, -PI/2);
        T0_p1p= T03p * T3_p1p;
        Tp1_4p = DHConvention(-PI/2, lu, 0, PI/2+t4); %
        T04p = T0_p1p * Tp1_4p;
        T45p = DHConvention(-PI/2, 0, 0, t5); 
        T05p = T04p * T45p;
        T5_p2p = DHConvention(-PI/2, 0, 0, PI/2); 
        T0_p2p = T05p * T5_p2p;
        Tp2_6p = DHConvention(0, lf, 0, t6); 
        T06p = T0_p2p * Tp2_6p;
        T67p = DHConvention(-PI/2, 0, 0, t7); 
        T07p = T06p * T67p;
        T78p= DHConvention(0, lw, 0, 0); 
        T08p=T07p*T78p;


T01= [ sin(t1), cos(t1), 0, 0
              0,               0, 1, 0
cos(t1), -sin(t1), 0, 0
              0,               0, 0, 1];

          
T02 =[sin(t1)*sin(t2), sin(t1)*cos(t2),  cos(t1), 0
       cos(t2),         -sin(t2),        0, 0
cos(t1)*sin(t2), cos(t1)*cos(t2), -sin(t1), 0
                    0,                       0,        0, 1];

                
                
T03= [cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2),   cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3), -cos(t2)*sin(t1), 0
                          cos(t2)*cos(t3),                            -cos(t2)*sin(t3),          sin(t2), 0
cos(t1)*cos(t3)*sin(t2) - sin(t1)*sin(t3), - cos(t3)*sin(t1) - cos(t1)*sin(t2)*sin(t3), -cos(t1)*cos(t2), 0
                                        0,                                           0,                0, 1];
                                    
                                    
                                   
 
                                    
                                    
                                    
T0_p1=[-cos(t2)*sin(t1), cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2),   cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3), 0
         sin(t2),                           cos(t2)*cos(t3),                            -cos(t2)*sin(t3), 0
-cos(t1)*cos(t2), cos(t1)*cos(t3)*sin(t2) - sin(t1)*sin(t3), - cos(t3)*sin(t1) - cos(t1)*sin(t2)*sin(t3), 0
               0,                                         0,                                           0, 1];  
  
           
           
T04=[cos(t2)*sin(t1)*sin(t4) - cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)), sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1), cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2), -lu*cos(t2)*sin(t1)
                                    cos(t2)*cos(t4)*sin(t3) - sin(t2)*sin(t4),                                   - cos(t4)*sin(t2) - cos(t2)*sin(t3)*sin(t4),                           cos(t2)*cos(t3),          lu*sin(t2)
cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4), cos(t1)*cos(t2)*cos(t4) - sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)), cos(t1)*cos(t3)*sin(t2) - sin(t1)*sin(t3), -lu*cos(t1)*cos(t2)
                                                                            0,                                                                             0,                                         0,                   1];           
 
                                                                        
                                                                        
T05= [- cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)), sin(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - cos(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)), sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1), -lu*cos(t2)*sin(t1)
                                                                - cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*sin(t5),                                                                 sin(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*cos(t5),                                   - cos(t4)*sin(t2) - cos(t2)*sin(t3)*sin(t4),          lu*sin(t2)
  cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)), cos(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)) - sin(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)), cos(t1)*cos(t2)*cos(t4) - sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)), -lu*cos(t1)*cos(t2)
                                                                                                                                              0,                                                                                                                                             0,                                                                             0,                   1];

                                                                                                                                          
                                                                                                                                          
T0_p2 = [- sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*cos(t4)*sin(t1),   cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)), sin(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - cos(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)), -lu*cos(t2)*sin(t1)
                                      cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4),                                                                   cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5),                                                                 sin(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*cos(t5),          lu*sin(t2)
  sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4), - cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) - sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)), cos(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)) - sin(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)), -lu*cos(t1)*cos(t2)
                                                                              0,                                                                                                                                               0,                                                                                                                                             0,                   1];

                                                                          
                                                                          
T06=[sin(t6)*(cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))) - cos(t6)*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)),   sin(t6)*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) + cos(t6)*(cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))), sin(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - cos(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)), - lf*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) - lu*cos(t2)*sin(t1)
                                                                                                    sin(t6)*(cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5)) + cos(t6)*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)),                                                                                                       cos(t6)*(cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5)) - sin(t6)*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)),                                                                 sin(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*cos(t5),                                               lf*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)) + lu*sin(t2)
cos(t6)*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - sin(t6)*(cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))), - sin(t6)*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - cos(t6)*(cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))), cos(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)) - sin(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)),   lf*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - lu*cos(t1)*cos(t2)
                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                   0,                                                                                                                                             0,                                                                                                         1];

                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                            
                                                                                                                                                                                                                                            
T07=[- cos(t7)*(cos(t6)*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) - sin(t6)*(cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)))) - sin(t7)*(sin(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - cos(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))), sin(t7)*(cos(t6)*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) - sin(t6)*(cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2)))) - cos(t7)*(sin(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) - cos(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))),   sin(t6)*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) + cos(t6)*(cos(t5)*(cos(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) - cos(t2)*sin(t1)*sin(t4)) + sin(t5)*(cos(t1)*sin(t3) + cos(t3)*sin(t1)*sin(t2))), - lf*(sin(t4)*(cos(t1)*cos(t3) - sin(t1)*sin(t2)*sin(t3)) + cos(t2)*cos(t4)*sin(t1)) - lu*cos(t2)*sin(t1)
                                                                                                                                                                      cos(t7)*(sin(t6)*(cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5)) + cos(t6)*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4))) - sin(t7)*(sin(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*cos(t5)),                                                                                                                                                                   - cos(t7)*(sin(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) - cos(t2)*cos(t3)*cos(t5)) - sin(t7)*(sin(t6)*(cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5)) + cos(t6)*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4))),                                                                                                       cos(t6)*(cos(t5)*(sin(t2)*sin(t4) - cos(t2)*cos(t4)*sin(t3)) + cos(t2)*cos(t3)*sin(t5)) - sin(t6)*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)),                                               lf*(cos(t4)*sin(t2) + cos(t2)*sin(t3)*sin(t4)) + lu*sin(t2)
  cos(t7)*(cos(t6)*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - sin(t6)*(cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)))) + sin(t7)*(sin(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) - cos(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))), cos(t7)*(sin(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) - cos(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))) - sin(t7)*(cos(t6)*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - sin(t6)*(cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2)))), - sin(t6)*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - cos(t6)*(cos(t5)*(cos(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) + cos(t1)*cos(t2)*sin(t4)) + sin(t5)*(sin(t1)*sin(t3) - cos(t1)*cos(t3)*sin(t2))),   lf*(sin(t4)*(cos(t3)*sin(t1) + cos(t1)*sin(t2)*sin(t3)) - cos(t1)*cos(t2)*cos(t4)) - lu*cos(t1)*cos(t2)
                                                                                                                                                                                                                                                                                                                                                                                                                      0,                                                                                                                                                                                                                                                                                                                                                                                                                     0,                                                                                                                                                                                                                                                   0,                                                                                                         1];
  

                                                                                                                                                                                                                                                                                                                                                                                                                  
%{                                                                                                                                                                                                                                                                                                                                                                                                                  
T1pp=simplify(expand(T01-T01p))
T2pp=simplify(expand(T02-T02p))
T3pp=simplify(expand(T03-T03p))
T4pp=simplify(expand(T04-T04p))
T5pp=simplify(expand(T05-T05p))
T6pp=simplify(expand(T06-T06p))   
T7pp=simplify(expand(T07-T07p))
Tp1pp=simplify(expand(T0_p1-T0_p1p))
Tp2pp=simplify(expand(T0_p2-T0_p2p))  
%}   
X1= [Lu;0;0;1];
Y1= simplify(T0_p1*X1)
X2= [Lf;0;0;1];
Y2= simplify(T0_p2*X2)