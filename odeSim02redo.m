function dxdt = odeSim02redo (t, input)
  dxdt = []; # fairly standard ODE45 script
  # k = [-0.081889,   0.051496,   0.041126,  -0.050160];
  # -0.081889
  #  0.051496
  #  0.041126
  # -0.050160

  # -0.080553
  #  0.050430
  #  0.040275
  # -0.049480

  k = [-0.080553, 0.050430, 0.040275, -0.049480];
  
  k11 = k(1); # could have skipped this way, not the most efficient
  k12 = k(2);
  k21 = k(3);
  k22 = k(4);
  dxdt(end+1) = k11*input(1)+k12*input(2);
  dxdt(end+1) = k21*input(1)+k22*input(2);

endfunction
