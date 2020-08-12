function mem = mem_der_LMN(v,i)
c = 5*10^(-7);
g_l = 2.5*10^(-8);
e_l = -0.07;
v_th = -0.05;
mem = (-1*(g_l/c)).*(v - e_l) + i./c;
end