function [v_final, spike_times,flag_out] = voltage_LMN(v_initial,i,N,flag_in)

e_l = -0.07;
v_th = -0.05;
v_r = -0.055;
h=0.001;
spike_times = zeros(N,1);
k1 = h*mem_der_LMN(v_initial,i);
k2 = h*mem_der_LMN(v_initial+k1, i);
v_final = v_initial +(k1 +k2)/2;
for n = 1:N
	if(flag_in(n,1) ~= 0)
		v_final(n,1) = v_r;
	end
end
flag_out = flag_in;
for n = 1:N
    if(v_final(n,1) > v_th)
        v_final(n,1) = v_r;
        spike_times(n,1) = 1;
        flag_out(n,1) = 1;
    end

%     if(v_final(n,1) < e_l)
%         v_final(n,1) = e_l;
%     end
    if(flag_in(n,1) == 1)
    	flag_out(n,1) = 2;
    end	
    if(flag_in(n,1) == 2)
    	flag_out(n,1) = 0;
    end	

end
