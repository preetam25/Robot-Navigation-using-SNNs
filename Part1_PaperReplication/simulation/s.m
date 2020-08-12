A = zeros(1,10000);
B = A;
A_c = A;
B_c = A;
A_f = 0;
B_f = 0;

d0 = 1;
motion = 1;
notmotion = 0;
A_c(1,1) = 3000e-12; %threshold current
B_c(1,1) = 2500e-12; %threshold current
A(1,1) = -0.07;
B(1,1) = -0.07;
w1 = 1e-13;
w2 = 1e-15;
w3 = 1e-15;
w4 = 1e-15;
spike_times_A = 0;
spike_times_B = 0;
for i=2:1:10000
    A_c(1,i) = 3000e-12 - w4*(A_f/i)*10000 + (B_f/i)*10000*w2;
    [A(1,i),spike_times_A] = voltage(A(1,i-1),A_c(1,i),1);
    A_f = A_f+spike_times_A;
    B_c(1,i) = 2500e-12 + w1*(A_f/i)*10000 - (B_f/i)*10000*w3;
    [B(1,i),spike_times_B] = voltage(B(1,i-1),B_c(1,i),1);
    B_f = B_f+spike_times_B;
end
plot(A)
hold on
plot(B)
hold off
    