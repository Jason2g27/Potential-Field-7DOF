(* ::Package:: *)

(* interactionProofMC.m *)
Clear[q1, q2, q3, q4, q5, q6, q7, lu, lf, lq, q1d, q2d, q3d, q4d, Ft, Fn]
(* DH transformation *)
DHC[alpha_, a_, d_, theta_] := {{Cos[theta], -Sin[theta], 0, a},
  {Sin[theta]*Cos[alpha], Cos[theta]*Cos[alpha], -Sin[alpha], -d*Sin[alpha]},
  {Sin[theta]*Sin[alpha], Cos[theta]*Sin[alpha], Cos[alpha], d*Cos[alpha]},
  {0, 0, 0, 1}};
  
(* Calculate Jacobian *)
CalculateJacobian[T01_, T02_, T03_, T04_, T05_, T06_, T07_, T08_ i_] := Module[
{z1, z2, z3, z4, z5, z6, z7, z8, z, pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8, pos, Jv},
  z1 = T01[[1 ;; 3, 3]];
  z2 = T02[[1 ;; 3, 3]];
  z3 = T03[[1 ;; 3, 3]];
  z4 = T04[[1 ;; 3, 3]];
  z5 = T05[[1 ;; 3, 3]];
  z6 = T03[[1 ;; 3, 3]];
  z7 = T04[[1 ;; 3, 3]];
  z8 = T05[[1 ;; 3, 3]];  
  z = {z1, z2, z3, z4, z5, z6, z7, z8};
  pos1 = T01[[1 ;; 3, 4]];
  pos2 = T02[[1 ;; 3, 4]];
  pos3 = T03[[1 ;; 3, 4]];
  pos4 = T04[[1 ;; 3, 4]];
  pos5 = T05[[1 ;; 3, 4]];
  pos6 = T06[[1 ;; 3, 4]];
  pos7 = T07[[1 ;; 3, 4]];
  pos8 = T08[[1 ;; 3, 4]];
  pos = {pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8};
  (* Initializing symbolic Jacobian matrix *)
  Jv = ConstantArray[0, {3, i - 1}];
  Do[
    Jv[[All, j]] = Cross[z[[j]], pos[[i]] - pos[[j]]],
    {j, i - 1}
  ];
  Return[Jv]
];

Assuming[
	(* Define variables *)
    Element[q1, Reals] && Element[q2, Reals] && Element[q3, Reals] && 
    Element[q4, Reals] && Element[q5, Reals] && Element[q6, Reals] && 
    Element[q7, Reals] && Element[lu, Reals] && Element[lf, Reals] && 
    Element[lw, Reals] && Element[q1d, Reals] && Element[q2d, Reals] && 
    Element[q3d, Reals] && Element[q4d, Reals] && Element[q5d, Reals] && 
    Element[q6d, Reals] && Element[q7d, Reals] && Element[Ft, Reals] && Element[Fn, Reals],
    (
        (* Homogeneous transformation *)
        %{
        T01 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(1,i)); %
        T12 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(2,i)); %
        T02 = T01*T12;
        T23 = DHConvention(pi/2, 0, 0, theta(3,i)); %
        T03 = T02*T23;
        T3_p1= DHConvention(-pi/2, 0, 0, -pi/2);
        T0_p1= T03 * T3_p1;
        Tp1_4 = DHConvention(-pi/2, lu, 0, pi/2+theta(4,i)); %
        T04 = T0_p1 * Tp1_4;
        T45 = DHConvention(-pi/2, 0, 0, theta(5,i)); 
        T05 = T04 * T45;
        T5_p2 = DHConvention(-pi/2, 0, 0, pi/2); 
        T0_p2 = T05 * T5_p2;
        Tp2_6 = DHConvention(0, lf, 0, theta(6,i)); 
        T06 = T0_p2 * Tp2_6;
        T67 = DHConvention(-pi/2, 0, 0, theta(7,i)); 
        T07 = T06 * T67;
        T78= DHConvention(0, lw, 0, 0); 
        T08=T07*T78;
        %}
        
		T01 = DHC[-\[Pi]/2, 0, 0, -\[Pi]/2 + q1];
		T12 = DHC[-\[Pi]/2, 0, 0, -\[Pi]/2 + q2];
		T23 = DHC[\[Pi]/2, 0, 0, q3];
		T3_p1 = DHC[-\[Pi]/2, 0, 0, -\[Pi]/2];
		Tp1_4 = DHC[-\[Pi]/2, lu, \[Pi]/2, q4];
		T34 = T3_p1 . Tp1_4;
		T45 = DHC[-\[Pi]/2, 0, 0, q5];
		T02 = T01 . T12;
		T03 = T02 . T23;
		T04 = T03 . T34;
		T05 = T04 . T45;
        
        T5_p2 = DHC[-\[Pi]/2, 0, 0, \[Pi]/2];
        Tp2_6 = DHC[0, lf, 0, q6];
        T56 = T5_p2 . Tp2_6;
        T67 = DHC[-\[Pi]/2, 0, 0, q7];
        T78 = DHC[0, lw, 0, q7];
        T06 = T05 . T56;
		T07 = T06 . T67;
		T08 = T07 . T78;

		Jve = CalculateJacobian[T01, T02, T03, T04, T05, T06, T07, T08 8];
		ve = Jve . {q1d, q2d, q3d, q4d, q5d, q6d, q7d};
		z8 = T08[[1 ;; 3, 3]]/Norm[T08[[1 ;; 3, 3]]];
		nt = Cross[ve, z4]/Norm(Cross[ve, z4]);
		
		Fnt = Ft * nt;
		torqt = Transpose[Jve] . Fnt;
		torqtsimp = Simplify[torqt];
		Fnn = Fn * z4;
		torqn = Transpose[Jve] . Fnn;
		torqnsimp = Simplify[torqn];
		Print["Tangential Torque: ", torqtsimp];
		Print["Normal Torque: ", torqnsimp];
		torqtmodified = Simplify[torqtsimp /. {q1 -> 0, q3 -> 0, q1d -> 0, q3d -> 0}];
		Print["Tangential Torque Sanity Test: ", torqtmodified]; 
		torqnmodified = Simplify[torqnsimp /. {q1 -> 0, q3 -> 0, q1d -> 0, q3d -> 0}];
		Print["Normal Torque Sanity Test: ", torqnmodified]; 
    )
]




