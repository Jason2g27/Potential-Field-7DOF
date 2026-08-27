
function Jv = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,T08,i)

    z1 = T01(1:3,3);
    z2 = T02(1:3,3);
    z3 = T03(1:3,3);
    z4 = T04(1:3,3);
    z5 = T05(1:3,3);
    z6 = T06(1:3,3);
    z7 = T07(1:3,3);
    z8 = T08(1:3,3);
    z = [z1,z2,z3,z4,z5,z6,z7,z8];
    pos1 = T01(1:3,4);
    pos2 = T02(1:3,4);
    pos3 = T03(1:3,4);
    pos4 = T04(1:3,4);
    pos5 = T05(1:3,4);
    pos6 = T06(1:3,4);
    pos7 = T07(1:3,4);
    pos8 = T08(1:3,4);
    pos = [pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8];
    %Jv = sym(zeros(3,i-1));
    Jv = zeros(3,7);
    for j = 1: i-1
        Jv(:,j) = cross(z(:,j), pos(:,i) - pos(:,j));
    end
end

