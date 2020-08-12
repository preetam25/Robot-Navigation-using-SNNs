function mem = mem_der_DTN(v,i)
c = 2*10^(-7);
g_l = 2.0*10^(-8);
e_l = -0.07;
v_th = -0.05;
mem = (-1*(g_l/c)).*(v - e_l) + i./c;
end