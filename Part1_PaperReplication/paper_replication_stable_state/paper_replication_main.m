tic
LMN = zeros(1000,2000);
DTN = zeros(1000,2000);
%mind that delta_t = 1ms
rotation =  -1*ones(1000,1); %sense of rotation
rotation(1:500,1) = -1*rotation(1:500,1);
% 1 means counterclockwise and -1 means clockwise
%first 500 -> ccw and last 500 -> cw in DTN
preferred_angle_LMN = zeros(1000,1);
preferred_angle_DTN = zeros(1000,1);
for i=1:1:1000
    preferred_angle_LMN(i,1) = 0.36*i;
    if(i<501)
        preferred_angle_DTN(i,1) = ((500 - i)/500) * 360;
    else
        preferred_angle_DTN(i,1) = ((i - 500)/500) * 360;
    end
end

% making the attractor network
weights_LMN_to_DTN = zeros(1000,1000); %excitatory signals
weights_DTN_to_LMN = zeros(1000,1000); %inhibitory weights
We_max = 63.3;
sigma_e = 67.2;
Wi_max = 59.5;
sigma_i = 147.3;

%w_kn = w_ij order of i and j with respect to paper for LMN_to_DTN
%w_nk = w_ij order of i and j with respect to paper for DTN_to_LMN
for i=1:1:1000
    for j=1:1:1000
        %delta1 = preferred_angle_LMN(j,1) - preferred_angle_DTN(i,1);
        if(i<501)
            delta1 = preferred_angle_LMN(j,1) - preferred_angle_DTN(i,1)+ 180 - 50;
            if(delta1 > 180)
                delta1 = delta1 - 360;
            elseif(delta1 < -180)
                delta1 = delta1 + 360;
            end
            % weights_LMN_to_DTN(i,j) = We_max * exp(-((delta1 + 180 - 50)^2)/(2*((sigma_e)^2)));
            weights_LMN_to_DTN(i,j) = We_max * exp(-((delta1 )^2)/(2*((sigma_e)^2)));
            
        elseif(i>500)
            delta1 = preferred_angle_LMN(j,1) - preferred_angle_DTN(i,1)+ 180 + 50;
            if(delta1 > 180)
                delta1 = delta1 - 360;
            elseif(delta1 < -180)
                delta1 = delta1 + 360;
            end
            %weights_LMN_to_DTN(i,j) = We_max * exp(-((delta1 + 180 + 50)^2)/(2*((sigma_e)^2)));
            weights_LMN_to_DTN(i,j) = We_max * exp(-((delta1 )^2)/(2*((sigma_e)^2)));
            
        end
        delta2 = preferred_angle_LMN(i,1) - preferred_angle_DTN(j,1);
        if(delta2 > 180)
            delta2 = delta2 - 360;
        elseif(delta2 < -180)
            delta2 = delta2 + 360;
        end
        
        weights_DTN_to_LMN(i,j) = Wi_max * exp(-((delta2)^2)/(2*((sigma_i)^2)));
        %now checking the threshold condition
        if(weights_LMN_to_DTN(i,j)<(0.55*We_max))
            weights_LMN_to_DTN(i,j) = 0;
        end
        if(weights_DTN_to_LMN(i,j)<(0.95*Wi_max))
            weights_DTN_to_LMN(i,j) = 0;
        end
    end
end

weights_DTN_to_LMN = -1*weights_DTN_to_LMN;

%making the integrator network
lateral_connection_weights = zeros(1000,500);
Wl_max = 25.6;
sigma_l = 35.1;
%w_li = w_ij order of i and j with respect to paper for lateral
for i=1:1:1000
    for j=1:1:500
        if(i<501)
            %delta = preferred_angle_DTN(j+500,1) - preferred_angle_DTN(i,1);
            delta = preferred_angle_DTN(j+500,1) - preferred_angle_DTN(i,1) - 30;
            
            if(delta > 180)
                delta = delta - 360;
            elseif(delta < -180)
                delta = delta + 360;
            end
            
            %lateral_connection_weights(i,j) = Wl_max * exp(-((delta - 30)^2)/(2*((sigma_l)^2)));
            lateral_connection_weights(i,j) = Wl_max * exp(-((delta )^2)/(2*((sigma_l)^2)));
            
        elseif(i>500)
            %delta = preferred_angle_DTN(501-j,1) - preferred_angle_DTN(i,1);
            delta = preferred_angle_DTN(501-j,1) - preferred_angle_DTN(i,1) + 30;
            
            if(delta > 180)
                delta = delta - 360;
            elseif(delta < -180)
                delta = delta + 360;
            end
            %lateral_connection_weights(i,j) = Wl_max * exp(-((delta + 30)^2)/(2*((sigma_l)^2)));
            lateral_connection_weights(i,j) = Wl_max * exp(-((delta )^2)/(2*((sigma_l)^2)));
            
        end
        if(lateral_connection_weights(i,j)<(0.60*Wl_max))
            lateral_connection_weights(i,j) = 0;
        end
    end
end

lateral_connection_weights = -1*lateral_connection_weights;

Fanin_LMN = cell(1000,1);
Fanin_DTN_DTN = cell(1000,1);
Fanin_DTN_LMN = cell(1000,1);

for i = 1:1000
    for j = 1:1000
        if(weights_DTN_to_LMN(i,j) ~= 0)
            Fanin_LMN{i} = [Fanin_LMN{i},j];
        end
    end
end


for i = 1:1000
    for j = 1:1000
        if(weights_LMN_to_DTN(i,j) ~= 0)
            Fanin_DTN_LMN{i} = [Fanin_DTN_LMN{i},j];
        end
    end
end

for i = 1:1000
    for j = 1:500
        if(lateral_connection_weights(i,j) ~= 0)
            if(i<501)
                Fanin_DTN_DTN{i} = [Fanin_DTN_DTN{i},j+500];
            else
                Fanin_DTN_DTN{i} = [Fanin_DTN_DTN{i},501-j];
            end
        end
    end
end

%define omega here
omega = 30*1e-9*ones(1,2000);

background_current_LMN = ones(1000,2000);
background_current_DTN = ones(1000,2000);

for i = 1:1000
    background_current_LMN(i,:) = rand*1200*(10^-9)*background_current_LMN(i,:);
end

for i = 1:1000
    background_current_DTN(i,:) = rand*(100*(10^-9)*background_current_DTN(i,:)+ 0.44*omega(1,:));
end

synaptic_current_LMN = zeros(1000,2000);
synaptic_current_DTN = zeros(1000,2000);

spike_times_LMN = zeros(1000,2000);
spike_times_DTN = zeros(1000,2000);
S = zeros(1000,2000);
for j = 1:10
    for i = 1:1000
        S(i,j) = 540*exp((-1*(i-500)^2)/(2*250*250));
    end
end
S = S*1e-9*0;

V_LMN = zeros(1000,2000);
V_DTN = zeros(1000,2000);
I_LMN = zeros(1000,2000);
I_DTN = zeros(1000,2000);

flag_LMN = zeros(1000,1);
flag_DTN = zeros(1000,1);
for time = 2:2000
    for LMN_i = 1:1000
        synaptic_current_LMN(LMN_i, time) = (1-0.01)*synaptic_current_LMN(LMN_i, time-1);
        for pre_neuron = 1:size(Fanin_LMN{LMN_i},2)
            synaptic_current_LMN(LMN_i, time) = synaptic_current_LMN(LMN_i, time) + 0.01*(10^-9)*spike_times_DTN(Fanin_LMN{LMN_i}(pre_neuron), time-1)*weights_DTN_to_LMN(LMN_i, Fanin_LMN{LMN_i}(pre_neuron));
        end
        synaptic_current_LMN(LMN_i, i) = synaptic_current_LMN(LMN_i, i) + 0.01*S(LMN_i,time);
    end
    
    for DTN_i = 1:1000
        synaptic_current_DTN(DTN_i, time) = (1-0.01)*synaptic_current_DTN(DTN_i, time-1);
        for pre_neuron = 1:size(Fanin_DTN_LMN{DTN_i},2)
            synaptic_current_DTN(DTN_i, time) = synaptic_current_DTN(DTN_i, time) + 0.01*(10^-9)*spike_times_LMN(Fanin_DTN_LMN{DTN_i}(pre_neuron), time-1)*weights_LMN_to_DTN(DTN_i, Fanin_DTN_LMN{DTN_i}(pre_neuron));
        end
        if(DTN_i < 501)
            for pre_neuron = 1:size(Fanin_DTN_DTN{DTN_i},2)
                synaptic_current_DTN(DTN_i, time) = synaptic_current_DTN(DTN_i, time) + 0.01*(10^-9)*spike_times_DTN(Fanin_DTN_DTN{DTN_i}(pre_neuron), time-1)*lateral_connection_weights(DTN_i, Fanin_DTN_DTN{DTN_i}(pre_neuron)-500);
            end
        else
            for pre_neuron = 1:size(Fanin_DTN_DTN{DTN_i},2)
                synaptic_current_DTN(DTN_i, time) = synaptic_current_DTN(DTN_i, time) + 0.01*(10^-9)*spike_times_DTN(Fanin_DTN_DTN{DTN_i}(pre_neuron), time-1)*lateral_connection_weights(DTN_i, 501 - Fanin_DTN_DTN{DTN_i}(pre_neuron));
            end
        end
    end
    %apply voltage and update spike
    I_LMN(:,time) = background_current_LMN(:,time) + synaptic_current_LMN(:,time);
    I_DTN(:,time) = background_current_DTN(:,time) + synaptic_current_DTN(:,time);
    [V_LMN(:,time), spike_times_LMN(:,time), flag_LMN] = voltage_LMN(V_LMN(:,time-1),I_LMN(:,time), 1000, flag_LMN);
    [V_DTN(:,time), spike_times_DTN(:,time), flag_DTN] = voltage_DTN(V_DTN(:,time-1),I_DTN(:,time), 1000, flag_DTN);
    
end

x = cell(1000,1);

for i = 1:1000
    for j = 1:2000
        if(spike_times_LMN(i,j)==1)
            x{i} = [x{i},j];
        end
    end
end

figure;
xlabel('Time in ms');
ylabel('Neuron ID');
for i = 1:1000
    scatter(x{i},i*ones(size(x{i})),3,'blue','fill');
    hold on;
end

toc



