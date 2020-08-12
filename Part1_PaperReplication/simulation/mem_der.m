function mem = mem_der(v,i)
c = 3*10^(-10);
g_l = 3*10^(-8);
e_l = -0.07;
v_th = 0.02;
mem = (-1*(g_l/c)).*(v - e_l) + i./c;
end