function [y, t] = step(obj,t, u, Ts)
[im_dot, iV_dot, iTwi, iTai] = deal(1,2,3,4); % Input indexes

mw_e_ = obj.mw_e;
dont_calculate_ABCD = 0;
x_ = obj.x;
Ts_c = mw_e_ / u(im_dot); 
assert(Ts_c <= obj.Ts_y_max, 'MMFVM: sampling period exceeded.');

switch obj.sampling_mode
    case 'equidistant'
        % Give output every (equidistant) sampling instance. Step model
        % when enough flow passed througt
        % --- Aggregate inputs
        u_sum_ = obj.u_sum + Ts * u;
        t_y_last_hit_ = obj.t_y_last_hit;
        
        y = obj.y_persist;        
        isOutputSampling = 0;
    case 'output'
        isOutputSampling = 1;
end

% --- When the aggegated volume exceeds element volume -> step MMFVM
while isOutputSampling || (u_sum_(im_dot) >= mw_e_)
    if isOutputSampling
        u_average = u;
        Ts_y = mw_e_ / u(im_dot);
        if t> 151
            1;
        end
    else 
        % --- Get "y" hit time and current "y" sampling time
        t_y_after_chunk = (u_sum_(im_dot) - mw_e_) / u(im_dot); %offset of "y sample" from current time
        t_y_hit = t + Ts - t_y_after_chunk;
        Ts_y = t_y_hit - t_y_last_hit_;
        
        % --- Get averaged inputs for the last Ts_y
        u_sum_after_chunk = t_y_after_chunk * u;
        u_sum_ = u_sum_ - u_sum_after_chunk;
        %   Aggregated water volume should now be exactly "mw_e"
        assert(u_sum_(im_dot) - mw_e_ < 1e-15);
        u_average = u_sum_ / Ts_y;
    end
    
    % --- Step the HX -----------------------------------------------------
    
    % --- Initialize state
    % x = [Twi, Tw, Two, Tb, Q]
    if isnan(x_(1))
        x_ = [...
            u(iTwi);...  % Twi;
            obj.Tw0;...  % Tw;
            0; ...       % Two
            obj.Tb0; ... % Tb
            0;];         % Q
    end
    
    % --- Get model matrices
    if ~dont_calculate_ABCD % no matrix change when more iterations in one step
        [A, B, C, D] = ABCD_exact_discretization(obj, u_average, Ts_y);
        %?! Precompute A, B?  Precompute their CF approximations?
        dont_calculate_ABCD = 1;
    end
    
    % --- Step state-space model (in "y" step)
    u_ = [u_average(iTwi); u_average(iTai)]; 
%     u_ = [u(iTwi); u(iTai)]; 
    
    %!!!! Enforce actual Twi in the state vector. Twi gets filled into
    %state vector with delay Ts_y!!!!!
    x_(1) = u_average(iTwi);

    x_ = A * x_ + B * u_; %~ x(k) = Ax(k-1) + Bu(k-1)!
    % all states delayed by Ts_y!!! due to Twi(k) getting in to x(1,k) = Twi(k)
    
    % outputing current value
    y = C * x_; % y(k) = Cx(k) + Du(k-1)????!! , but D are zeros!
    
    if isOutputSampling
        break; % Do only one run through the while loop
    else
        % --- Actualize sampling state
        t_y_last_hit_ = t_y_hit;
        u_sum_ = u_sum_after_chunk;
    end
end

% --- Step simulation time (in "u" steps)
if isOutputSampling
    t = t + Ts_y;
else
    t = t + Ts;
end

% --- Save state
obj.x = x_;
if ~isOutputSampling
    obj.y_persist = y;
    obj.u_sum = u_sum_;
    obj.t_y_last_hit = t_y_last_hit_;
end
end