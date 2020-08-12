clear all;
tic;

%% Defining constants

Rp = 2e-3; % refractory period in s
w = 15e3;
wa = w;
wm = w;
ws = w;
wo = w;
lambda = 100; % spike_times/s
I0 = 1e-12; % in A
dt = 1e-3; % time step in s
T = 0.5; % total simulation time in s
VT = 0.020; % threshold voltage in V
EL = -0.070; % resting potential in V
tm = 0.015; % time constant of membrane in s
ts = tm/4; % synaptic time constant in s
N = 9;
dc = 1e5*I0;
VTL = -0.1;
window = 10;

%% Bot control
% neuron 1 - position A
% neuron 2 - position B
% neuron 3 - position C
% neuron 4 - position D
% neuron 5 - angle 0
% neuron 6 - angle 90
% neuron 7 - angle 180
% neuron 8 - angle 270
% neuron 9 - motor

potentials = EL*ones(N,T/dt);
currents = zeros(N,T/dt);
spikes = zeros(N, T/dt);

t1 = 50; %A to B begin
t2 = 60; %A to B end

t3 = 150; %B to C begin
t4 = 160; %B to C end

t5 = 250; %C to D begin
t6 = 260; %C to D end

t7 = 350; %D to A begin
t8 = 360; %D to A end

t9 = 450; %A to B begin
t10 = 460; %A to B end

currents(9,t1:t2) = dc;
currents(9,t3:t4) = dc;
currents(9,t5:t6) = dc;
currents(9,t7:t8) = dc;
currents(9,t9:t10) = dc;

r1 = randi([1,50],1,1);
r2 = randi([1,50],1,1);
r3 = randi([1,50],1,1);
r4 = randi([1,50],1,1);
r5 = randi([1,50],1,1);

currents(5,1:t2+r1) = dc;
currents(6,t2+r1+1:t4+r2) = dc;
currents(7,t4+r2+1:t6+r3) = dc;
currents(8,t6+r3+1:t8+r4) = dc;
currents(5,t8+r4+1:T/dt) = dc;


%% Solving the network

for t = 1:T/dt-1
        
    if(t<=window)
        currents(1,t) = dc;
    else
        currents(1,t) = 0;
    end
    
    currents(2,t) = 0;
    currents(3,t) = 0;
    currents(4,t) = 0;
    
    if(t>window)
        for z = t-window:t
            element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));

            if(spikes(1,z)==1)
                currents(1,t) = currents(1,t) + ws*element;
            end

            if(spikes(2,z)==1)
                currents(2,t) = currents(2,t) + ws*element;
            end

            if(spikes(3,z)==1)
                currents(3,t) = currents(3,t) + ws*element;
            end

            if(spikes(4,z)==1)
                currents(4,t) = currents(4,t) + ws*element;
            end
        end
        
       if(sum(spikes(9,t-window:t),'all')~=0)
           
            if(sum(spikes(5,t-window:t),'all')~=0)
                if(sum(spikes(1,t-window:t),'all')~=0)
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(5,z)==1)
                            currents(2,t) = currents(2,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(2,t) = currents(2,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(2,t) = currents(2,t)+wo*element;
                        end
                        currents(1,t) = - currents(2,t);
                    end
                elseif(sum(spikes(4,t-window:t),'all')~=0)    
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(5,z)==1)
                            currents(3,t) = currents(3,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(3,t) = currents(3,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(3,t) = currents(3,t)+wo*element;
                        end
                        currents(4,t) = - currents(3,t);
                    end
                end
            end
            
            if(sum(spikes(6,t-window:t),'all')~=0)
                if(sum(spikes(2,t-window:t),'all')~=0)
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(6,z)==1)
                            currents(3,t) = currents(3,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(3,t) = currents(3,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(3,t) = currents(3,t)+wo*element;
                        end
                        currents(2,t) = - currents(3,t);
                    end
                elseif(sum(spikes(1,t-window:t),'all')~=0)    
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(6,z)==1)
                            currents(4,t) = currents(4,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(4,t) = currents(4,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(4,t) = currents(4,t)+wo*element;
                        end
                        currents(1,t) = - currents(4,t);
                    end
                end
            end
            
            if(sum(spikes(7,t-window:t),'all')~=0)
                if(sum(spikes(2,t-window:t),'all')~=0)
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(7,z)==1)
                            currents(1,t) = currents(1,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(1,t) = currents(1,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(1,t) = currents(1,t)+wo*element;
                        end
                        currents(2,t) = - currents(1,t);
                    end
                elseif(sum(spikes(3,t-window:t),'all')~=0)    
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(7,z)==1)
                            currents(4,t) = currents(4,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(4,t) = currents(4,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(4,t) = currents(4,t)+wo*element;
                        end
                        currents(3,t) = - currents(4,t);
                    end
                end
            end
            
            if(sum(spikes(8,t-window:t),'all')~=0)
                if(sum(spikes(3,t-window:t),'all')~=0)
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(8,z)==1)
                            currents(2,t) = currents(2,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(2,t) = currents(2,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(2,t) = currents(2,t)+wo*element;
                        end
                        currents(3,t) = - currents(2,t);
                    end
                elseif(sum(spikes(4,t-window:t),'all')~=0)    
                    for z = t-window:t
                        element = I0*(exp(-(t-z)*dt/tm) - exp(-(t-z)*dt/ts));
                        if(spikes(8,z)==1)
                            currents(1,t) = currents(1,t)+wa*element;
                        end
                        if(spikes(9,z)==1)
                            currents(1,t) = currents(1,t)+wm*element;
                        end
                        if(spikes(1,z)==1)
                            currents(1,t) = currents(1,t)+wo*element;
                        end
                        currents(4,t) = - currents(1,t);
                    end
                end
            end
        end 
    end
    
    if(t>=3)
        [potentials(:,t-1:t+1), spikes(:,t)] = LIF(currents(:,t-1), currents(:,t), potentials(:,t-2:t), t, dt, VT, EL, VTL);
    end
end

%% Plotting firing pattern

t = 1:T/dt;
x = cell(N,1);
 
 for i = 1:N
    for j = 1:T/dt
        if(spikes(i,j)==1) 
            x{i} = [x{i},j]; 
        end
    end
 end
 
 figure;
 colors = ['r','g','b','y','k','k','k','k','m'];
 for i = 1:N
    scatter(x{i},i*ones(size(x{i})),10,colors(i),'fill');
    hold on;
 end
 xlabel('Time in ms'); 
 ylabel('Neuron ID'); 
 
figure;
plot(t, potentials(1,:), 'r', t, potentials(2,:), 'g', t, potentials(3,:), 'b', t, potentials(4,:), 'y');
title('Time evolution of membrane potentials');
xlabel('Time (ms)');
ylabel('Membrane potential (V)');
legend('neuron A','neuron B','neuron C','neuron D');

toc;

% add lower threshold