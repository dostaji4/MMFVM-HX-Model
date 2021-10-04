function y = ss(obj,u)
    
    UAwb = obj.UAwb_model(m_dot*3600);
    UAba = obj.UAba_model(V_dot*3600);

    % --- Get model matrices ---
    [A,B,C,D] = obj.ABCD_fvm_discrete(m_dot,...
        obj.Nw, obj.Nwb, obj.mw, obj.Cw, obj.Cb, UAwb, UAba, obj.Ts_max);

    % --- Compute steady state ---
    Tss = (eye(size(A)) - A)\B*[Twi; Tai];
    yss = C*Tss + D*[Twi; Tai];

    Twss = Tss(1:obj.Nw);
    Tbss = Tss(obj.Nw+1:end);

    if nargout == 0
        % Plot
        x_edge = linspace(0,1,obj.Nw+1);
        x_mid  = x_edge(1:end-1)+ diff(x_edge(1:2))/2;
        p = plot(x_mid, Tbss,'+-',...
            x_mid, Twss,'+-',...
            'LineWidth',1);
        p(1).DisplayName = 'B fvm_d';
        p(2).DisplayName = 'W fvm_d';
        legend('off');
        legend('show');
    end
end