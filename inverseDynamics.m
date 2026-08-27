function tau = inverseDynamics(state, M,C,G)
    tau = M * state(15:21) + C + G;

end