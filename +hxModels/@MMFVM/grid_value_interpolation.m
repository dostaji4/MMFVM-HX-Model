function [Tw_output, Two_output] = grid_value_interpolation(Tw,p,schemeIdx,varargin)
% Interpolate state values at grid points from moving water elements
% SchemeIdx
%   0 - constant upwind (first order accurate)
%   1 - linear (second order)
%   2 - quadratic (third order accurate)
%   3 - TVD (linear + quadratic + flux limiter)

isPlot = 0;
if nargin > 3
    z = varargin{1}(:);
    %     z = z(:);
    Nw = numel(z);
end

if nargin > 4 % plot
    t = varargin{2};
    plotTime_start = inf;
    plotTime_end = 100;
    isPlot =  plotTime_start < t & t < plotTime_end;
end

giveApprox = 0;
if nargin > 5 % give approximation at zn
    zn = varargin{3}(:);
    giveApprox = 1;
end

Two_position_exact = 0; % 1 ~ Two = T('x=1')
%% Out
dz = z(2)-z(1);
if numel(Tw) == numel(z) + 2 %MMFVM input    
    z_static = [-dz/2;z];
    z_actual = z_static + p*dz;
    Two = Tw(end);
    Tw = Tw(1:end-1);
elseif numel(Tw) == numel(z) %FVM input
    z_static = z;
    z_actual = z;
    Two = Tw(end);
else
    error('Invalid input.');
end
   
switch schemeIdx
    case 0 % No interpolation (First order upwind)
        [~,f_Tw] = const_appr(z_actual,Tw,z);
    case 1 % Linear interpolation (Second order central)
        [~,f_Tw] = linear_appr(z_actual,Tw,z);
    case 2 % Quadratic interpolation (Third order upwind)
        
        [~,f_phi3u_z] = phi_qm(z_actual,Tw,z,0); %upwind
        [~,f_phi3d_z] = phi_qm(z_actual,Tw,z,1); %downwind
        f_phi3_cmb_z = @(z) max([f_phi3u_z(z),f_phi3d_z(z)],[],2);
        
        [~,f_Tw2_z] = linear_appr(z_actual,Tw,z);
        [~,f_Tw3u_z] = quadr_appr(z_actual,Tw,z,0);%upwind
        [~,f_Tw3d_z] = quadr_appr(z_actual,Tw,z,1);%downwind
        
        f_Tw3_z = @(z) (1-f_phi3_cmb_z(z)).*f_Tw2_z(z)...
            + f_phi3_cmb_z(z).*(f_phi3u_z(z).*f_Tw3u_z(z) +...
            f_phi3d_z(z).*f_Tw3d_z(z))./...
            (f_phi3u_z(z)+f_phi3d_z(z)+eps);
        
        f_Tw = f_Tw3_z;
    case 3 % TVD
        %{
        !! TODO f_Tw(z)
        % --- Smoothness monitor
        % Jakobsen, Chemical Reactor Modelling, p. 1140
        % r(i) = (t-T(i-1)) / (T(i+1)-t)
        r = (T_P - T_W) ./ (T_E - T_P); % forward monitor
        %         r = (T_P - T_E) ./ (T_W - T_P); % backwards monitor
        
        
        % --- Flux limiter
        phi = ( r+abs(r)) ./ (1+r); %VanLeer
        %                                 phi = max(0,min(1,r,'includenan'),'includenan'); %minmod
        %                                 phi = max(0,min(2,min(2*r,(1+r)/2),'includenan'),'includenan'); %MUCSL
        %                                 beta = 2; % beta in [1,2]
        %                                 phi = max(0, min(beta, min(beta*r, (1+r)/2,'includenan'),'includenan'),'includenan'); %Gen. minmod
        %                                 phi = max(0, max(min(beta*r,1),min(r, beta))); %Sweby
        %                                  phi = max(0,2*r ./ (r.^2 + 1),'includenan'); %albada 2
        %                                 phi = max(0,max(min(2*r,1,'includenan'), min(r,2,'includenan'),'includenan'),'includenan'); %superbee
        %                                 beta = 100; % width of window in (2,inf)
        %                                 phi = max(0,min(1,beta/2-abs(r-beta/2),'includenan'),'includenan'); %dost2
        
        
        phi = [0; phi]; %add phi for Tw(1)
        % Disable ill conditioned H approximation
        phi = phi .* phi_Hd;
        phi(isnan(phi)) = 0;
        
        f_TwTwo_3 = @(p) (1-repmat(phi,1,numel(p))).*f_1u(p) +...
            repmat(phi,1,numel(p)).*f_TwTwo2c(p);
        f_TwTwo = @(p) f_TwTwo_3(p);
        if 0 &&  t > 50
            
            pp = linspace(0,1,100);
            dz = 1/Nw;
            zp = z;
            yp = f_TwTwo_3(fliplr(pp));
            xp = zp + dz.*pp;
            figure(9);
            plot(xp',yp(1:end-1,:)');
            pause(0.05);
        end
        
        TwTwo_2c = f_TwTwo2c(p);
        TwTwo_3u = [TwTwo_2c(1);f_3u(p)];
        TwTwo_3d = [f_3d(p);TwTwo_2c(end)];
        TwTwo_TVD = (1-phi) .* TwTwo_2c + ...
            phi .* TwTwo_3u;
        
        Tw = TwTwo_TVD(1:end-1);
        
        % SPECIAL CARE FOR TWO (essectially H, but
        % when there are bumps, then L)
        Two = TwTwo_2c(end);
        
        % Plot
        if isPlot
            plot_3tvd;
        end
        %}
end

Tw_output = f_Tw(z);
if Two_position_exact
    Two_output = f_Tw(1);    
else
    Two_output = f_Tw(1+dz/2);
    Two_output = p*Tw(end) + (1-p)*Two; % interpolate linearly to get end value
end

if giveApprox
    Tw_output = f_Tw(zn);
    
    % --- Modify the trailing end for smooth transition
    idx_e = zn > z_static(end); % select the trailing end 
    ne = nnz(idx_e);
    ze_old = zn(idx_e); % old positions 
    ze = linspace(ze_old(1),z_actual(end)+dz,ne); % new positions 
    Two_ = p*Tw(end) + (1-p)*Two; % interpolate linearly to get end value
    Twe = interp1([ze(1);z_actual(end)+dz],... %from first of the trailing end to the actual position of last old element
        [f_Tw(ze(1));Two_],... %from current interpolation into interpolated end value at its position
        ze,'linear'); % interpolate linearly on the new positions
    Tw_output(idx_e) = Twe; %replace the end
    
%{
    figure(7);cla    
    plot_profile(f_Tw);    
    plot(zn,Tw_output,'+');
    plot(z_actual(end)+dz,Two,'*');
    plot(ze,Twe,'.');
%     pause
%}
end

if isPlot
    figure(7);cla    
    plot_profile(f_Tw);
    pause(0.1);
end

    function plot_profile(h_f_Tw)       
       plot_profile_base(z,p,20,70,1);
       x = linspace(0,z_static(end)+dz,1000);
       y = h_f_Tw(x);
       plot(x,y);       
       plot(z_actual,Tw,'*'); hold on    
    end

%% Constant approximation (First order upwind)
    function [yn, h_f_yn] = const_appr(z,y,zn)
        yn = f_yn(zn);
        h_f_yn = @f_yn;
        function yn = f_yn(zn)
            % behavior on faces - upwind                
            dz = diff(z);
            assert(all((dz - dz(1))<eps));
            dz = dz(1);
            z_ext = [z(:)+dz/2];        
            assert(numel(z) == numel(y));
            z_ext_ = repmat(z_ext(:),1,numel(zn));
            zn_ext_ = repmat(zn(:)',numel(z),1);
            d = zn_ext_ <= z_ext_;
            d = xor([d; zeros(1,size(d,2))],[zeros(1,size(d,2));d]);
            d = d(1:end-1,:);
            yn = d'*y;
        end
    end
%% Linear approximation (Second order central - low precision, high resolution flux)
% Two = T(x=x_Two)
% ~ artificial delay in Two

    function [yn,h_f_yn] = linear_appr(z,y,zn)
        f_yn = @(zn) interp1(z(:),y(:),zn(:),'linear','extrap');
        yn = f_yn(zn);
        h_f_yn = f_yn;
    end

%% 3rd order (high precision, low resolution flux)
    function [yn,h_f_yn,xh] = quadr_appr(z,y,zn, dir)
        
        if nargin < 4
            dir = 0; % UPWIND
        end
        assert(any(dir==[0,1]),'Unknown direction setting. Give 0 for upwind (default) or 1 for downwind.');
        assert(numel(z) == numel(y),'"z" and "y" do not have the same dimension.');
        assert(numel(y) >= 3,'Quadratic approximation needs at least 3 original points.');
        
        yW = y(1:end-2);
        yP = y(2:end-1);
        yE = y(3:end);
        a = (yE-2*yP+yW)/2;
        b = (yE-yW)/2;
        c = yP;
        
        yn = f_yn(zn);
        h_f_yn = @f_yn;
        xh = -(b./(2*a));  %parabola head relative position (to P, in units dz)
        
        function yn = f_yn(zn)            
            f_yn_p = @(p,ind) a(ind)*p.^2 + b(ind)*p + c(ind);
            dz = diff(z);
            assert(all(dz - dz(1) < eps),'Only valid for uniform grids.');
            dz = dz(1);
            yn = NaN(numel(zn),1);
            % Get all
            for i = 1:numel(zn)
                if zn(i) < z(1)
                    % extrapolate behind W point
                    pzn = (z(1)-zn(i))/dz; %relative position to W point
                    yn(i) = f_yn_p(-1-pzn,1); % p=-1~W,p=0~P,p=1~E
                elseif zn(i) > z(end)
                    % extrapolate behind E point
                    pzn = (zn(i)-z(end))/dz; %relative position from E point
                    yn(i) = f_yn_p(1+pzn,numel(a)); % p=-1~W,p=0~P,p=1~E
                else
                    if ~dir %upwind
                        if zn(i) <= z(2)
                            % interpolate between W and P points (NOT UPWIND HERE)
                            pzn = (z(2)-zn(i))/dz; %relative position to P (nearest greater)
                            yn(i) = f_yn_p(-pzn,1);% p=-1~W,p=0~P,p=1~E
                        else
                            % interpolate between P and E points
                            i_ng = find(zn(i) <= z,1); %nearest greater (E)
                            pzn = (z(i_ng)-zn(i))/dz; %relative position to E point
                            yn(i) = f_yn_p(1-pzn,i_ng-2);% p=-1~W,p=0~P,p=1~E
                        end
                    else %downwind
                        if zn(i) >= z(end-1)
                            % interpolate between P and E points (NOT DOWNWIND HERE)
                            pzn = (zn(i)-z(end-1))/dz; %relative position from P (nearest smaller)
                            yn(i) = f_yn_p(pzn,numel(a));% p=-1~W,p=0~P,p=1~E
                        else
                            % interpolate between W and P points
                            i_ng = find(zn(i) < z,1); %nearest greater (P)
                            pzn = (z(i_ng)-zn(i))/dz; %relative position to P point
                            yn(i) = f_yn_p(-pzn,i_ng-1);% p=-1~W,p=0~P,p=1~E
                        end
                    end
                    
                end
            end
        end
    end

    function r = rate_monitor(y)
        r = -1*ones(numel(y),1); % default value for r(1) and r(end). -1 gives parabola head = 0
        yW = y(1:end-2);
        yP = y(2:end-1);
        yE = y(3:end);
        r(2:end-1) = (yP-yW) ./ (yE - yP);
    end

    function [phi,h_phi3_z] = phi_qm(z,y,zn,dir)        
        assert(numel(y) == numel(z));        
        if nargin < 4
            dir = 0; % UPWIND
        end
        % soft switch
        x_switch = 1.1;
        beta = 1;
        f_xh_r = @(r) (r+1)./(2*(r-1)); % parabola head position given by r = (P-W)/(E-P)
        logistic = @(xh) beta./(1+exp(-7*(xh-x_switch)));
        f_phi3_r = @(r) logistic(abs(f_xh_r(r)))-logistic(0);
        
        f_phi3_y = @(y) f_phi3_r(rate_monitor(y));
        
        phi = f_phi3_z(zn,dir);
        h_phi3_z = @(z)f_phi3_z(z,dir);
        
        %extend for zn
        function phi = f_phi3_z(zn,dir)
            phi3_y = f_phi3_y(y);
            phi = zeros(numel(zn),1);
            for i = 1:numel(zn)
                if ~dir %Upwind
                    idx = find(zn(i) > z,1,'last');
                    if ~isempty(idx)
                        if idx ~= numel(phi3_y)
                            phi(i) = phi3_y(idx);
                        else
                            phi(i) = phi3_y(end-1); % HACK! allow for extrapolation
                        end
                    end
                else % Downwind
                    idx = find(zn(i) < z,1);
                    if ~isempty(idx)
                        phi(i) = phi3_y(idx);
                    end
                end
            end
        end
    end

%% ============================== PLOT ====================================

    function plot_2uw
        
        Tmax = 80;Tmin = 20;
        zp = plot_profile_base(z,p,Tmin,Tmax,Two_position_exact);
        
        % --- Plot current element temperatures
        Tw_ = TwTwo(1:end-1);
        plot(zp+dz*p,Tw_,'x');
        
        % --- Plot interpolation functions
        xp = linspace(0,1,2);
        yp = f_TwTwo2c(xp);
        xp_ = repmat(zp(2:end),1,numel(xp)) + dz*(p +...
            repmat(-xp,Nw,1));
        plot(xp_',yp(1:end-1,:)','LineWidth',1);
        
        % --- Plot Tw output value
        TwTwo_2c = f_TwTwo2c(p);
        plot(zp(2:end),TwTwo_2c(1:end-1),'o');
        
        % --- Plot Two value
        if Two_position_exact
            plot(1,f_Two2(p),'o');
            if zp(end)+p*dz < 1
                xp = linspace(p,0.5,10);
                yp = f_Two2(xp);
                plot(1+(p-xp)*dz,yp,'--','Linewidth',1);
            end
            %             plot(1,f_Two2u(p),'o');
            %             xp = linspace(p,1+0.5,10);
            %             yp = f_Two2u(xp);
            %             plot(1+(p-xp)*dz,yp,'--','LineWidth',1)
        else
            %plot(zp(end),f_Two2(p),'o');
            % nothing to plot
        end
    end

    function plot_3uw
        Tmax = 80;Tmin = 20;
        zp = plot_profile_base(z,p,Tmin,Tmax,Two_position_exact);
        
        % --- Plot current element temperatures
        plot(zp+dz*p,TwTwo,'x');
        
        % --- Plot interpolation
        fp = @(x) a*(1-x).^2 + b*(1-x) + repmat(c,1,numel(x));
        xp = linspace(0,2,25);
        yp = fp(xp);
        xp_ = repmat(zp(3:end),1,numel(xp)) + dz*(p +...
            repmat(-xp,Nw,1));
        clear h
        h = plot(xp_(1:end-1,:)',yp(1:end-1,:)','LineWidth',1);
        
        % --- Plot output fcn
        TwTwo_3u = [TwTwo_2c(1);f_3u(p)];
        Twp = (1-phi_Hu).*TwTwo_2c + (phi_Hu).*TwTwo_3u;
        clear l
        if Two_position_exact
            zp(end) = 1;
        end
        for j = 1:numel(Twp)
            
            l = plot(zp(j+1),Twp(j),'o');
            if j >=2 && j <= numel(h)+1
                l.Color = h(j-1).Color;
            end
            markersize = 8;
            if phi_Hu(j) < 0.3
                l.Marker = '+';
                l.MarkerSize = markersize;
                l.Color = 'k';
            else
                l.MarkerSize = phi_Hu(j)*markersize;
            end
            
        end
        title('2Upwind');
        drawnow;
    end

    function plot_3dw
        Tmax = 80;Tmin = 20;
        zp = plot_profile_base(z,p,Tmin,Tmax,Two_position_exact);
        
        % --- Plot current element temperatures
        plot(zp+dz*p,TwTwo,'x');
        
        % --- Plot interpolation
        fp = @(x) a*(x-1).^2 + b*(x-1) + repmat(c,1,numel(x));
        xp = linspace(0,2,25);
        yp = fp(xp);
        xp_ = repmat(zp(3:end),1,numel(xp)) + dz*(p +...
            repmat(-xp,Nw,1));
        clear h
        h = plot(xp_(1:end-1,:)',yp(1:end-1,:)','LineWidth',1);
        
        % --- Plot output fcn
        TwTwo_2c = f_TwTwo2c(p);
        TwTwo_3d = [f_3d(p);TwTwo_2c(end)];
        Twp = (1-phi_Hd).*TwTwo_2c + (phi_Hd).*TwTwo_3d;
        clear l
        if Two_position_exact
            zp(end) = 1;
        end
        for j = 1:numel(Twp)
            
            l = plot(zp(j+1),Twp(j),'o');
            if j <= numel(h)
                l.Color = h(j).Color;
            end
            markersize = 8;
            if phi_Hd(j) < 0.3
                l.Marker = '+';
                l.MarkerSize = markersize;
                l.Color = 'k';
            else
                l.MarkerSize = phi_Hd(j)*markersize;
            end
            
        end
        title('2Downwind');
        drawnow;
    end

    function plot_3uwdw
        
        Tmax = 80;Tmin = 20;
        zp = plot_profile_base(z,p,Tmin,Tmax,Two_position_exact);
        
        % --- Plot current element temperatures
        plot(zp+dz*p,TwTwo,'x');
        
        % --- Plot interpolation
        fp = @(x) a*(x-1).^2 + b*(x-1) + repmat(c,1,numel(x));
        xp = linspace(0,2,25);
        yp = fp(xp);
        xp_ = repmat(zp(3:end),1,numel(xp)) + dz*(p +...
            repmat(-xp,Nw,1));
        clear h
        h = plot(xp_(1:end-1,:)',yp(1:end-1,:)','LineWidth',1);
        
        % --- Plot output fcn
        TwTwo_2c = f_TwTwo2c(p);
        TwTwo_3u = [TwTwo_2c(1);f_3u(p)];
        TwTwo_3d = [f_3d(p);TwTwo_2c(end)];
        phi_combined = max(phi_Hd,phi_Hu);
        Twp = (1-phi_combined).*TwTwo_2c + phi_combined .* ...
            (phi_Hd.*TwTwo_3d + phi_Hu.*TwTwo_3u)./...
            (phi_Hd+phi_Hu);
        clear l
        if Two_position_exact
            zp(end) = 1;
        end
        for j = 1:numel(Twp)
            
            l = plot(zp(j+1),Twp(j),'o');
            if j <= numel(h)
                l.Color = h(j).Color;
            end
            markersize = 8;
            if phi_combined(j) < 0.3
                l.Marker = '+';
                l.MarkerSize = markersize;
                l.Color = 'k';
            else
                l.MarkerSize = phi_combined(j)*markersize;
            end
            
        end
        drawnow;
    end

    function plot_3tvd
        figure(7);
        yp = [T_W(end),T_P(end),T_E(end)];
        xp = [-1, 0, 1] -0.5 + p;
        
        % full line
        xp_ = linspace(xp(1)-1,xp(3)-1,15);
        fp = @(xp_) (T_E(end)-2*T_P(end)+T_W(end))/2 * (1.5-(p-xp_)).^2 + ...
            (T_E(end)-T_W(end))/2 * (1.5-(p-xp_)) + ...
            T_P(end);
        yp_ = fp(xp_);
        pl = plot(xp_+1,yp_,'LineWidth',1,'LineStyle','-');
        ax = gca;
        ax.XTick = [-1.5,-1,-0.5,0,0.5,1];
        ax.XTickLabel = {'W','','P','','E','1'};
        hold on
        % extrapolate
        if xp(3)-1 < 0
            xp_ = linspace(xp(3)-1,0,10);
            yp_ = fp(xp_);
            plot(xp_+1,yp_,'LineWidth',1,'LineStyle','--','Color',pl.Color);
        end
        % Nodes
        TwTwo_2c = f_TwTwo2c(p);
        plot(xp,yp,'x','Color',pl.Color);
        % Second order
        plot(1,TwTwo_2c(end),'+','Color',pl.Color);
        % Trird order approx.
        plot(1,TwTwo_H(end),'x','Color',pl.Color);
        % Final
        plot(1,Two,'o','Color',pl.Color);
        % Lines
        %                                     plot([xp; xp]+0.5-p,repmat([min(yp);max(yp)],1,3),'k--','LineWidth',1);
        
        %                                     %Plot boxes
        %                                     boxx = [0 1 1 0 0]-0.5; boxy = [0 0 1 1 0]-0.5;
        %                                     for j = 1:numel(xp)
        %                                         patch('XData',xp(j)+boxx,'YData',yp(j)+boxy,'FaceAlpha',0.1*p,'LineWidth',1,'EdgeColor',[1 1 1]*0.9-0.1*p);
        %                                     end
        
        h = get(figure(7),'UserData');
        if ~isempty(h), delete(h);end
        h = annotation('textbox',[.05 .0 .0 .95],'String',{['p=' num2str(p)];['t=' num2str(t)]},'FitBoxToText','on');
        set(figure(7),'UserData',h);
        drawnow;
    end

    function zb = plot_profile_base(z,p,ymin,ymax,Two_exact)
        sy = ymax-ymin;
        oy = ymin;
        box_height = sy/10;
        n = numel(z);
        dz = z(2) - z(1);
        
        % --- Plot water lines
        z_ = [-dz/2; z];
        cl = 0.8*[0.8 0.8 1];
        yl_min = ymin;
        yl_max = ymax;
        plot([z_,z_]'+p*dz,...
            [yl_min * ones(n+1,1), yl_max * ones(n+1,1)]',...
            'LineStyle','--','Color',cl,'LineWidth',1); hold on
        
        % ---  Plot body lines
        if Two_exact
            z_ = [z; 1];
        else
            z_ = [z; 1+dz/2];
        end
        cl = 0.8*[1 1 1];
        yl_min = ymin - box_height;
        yl_max = ymax;
        plot([z_,z_]',...
            [yl_min * ones(n+1,1), yl_max * ones(n+1,1)]',...
            'LineStyle','--','Color',cl,'LineWidth',1); hold on
        
        
        ax = gca;
        ax.XGrid = 'off';
        
        zb = [-dz/2; z]; %add Twi
        box_height = sy/10;
        plot_boxes(zb,p,box_height,oy);
    end

    function plot_boxes(zp,p,y_height,y_offset)
        % --- Plot boxes
        boxx = [0 1 1 0 0]-0.5; boxy = [0 0 1 1 0]-0.5;
        y0w = y_offset;
        y0b = y_offset - y_height;
        cb = [1 1 1]*0.6;
        cw = [0.8 0.8 1]*0.8;
        dz = median(diff(zp));
        % Tb
        for j = 2:numel(zp)
            
            patch('XData',zp(j)+(boxx)*dz,'YData',y0b+boxy*y_height,...
                'FaceColor',cb,'FaceAlpha',0.1,...
                'LineWidth',1,'EdgeColor',cb);
            plot(zp(j)+(0)*dz,y0b,'Marker','+','Color',cb);
        end
        
        %Tw
        for j = 1:numel(zp)
            patch('XData',zp(j)+(boxx+p)*dz,'YData',y0w+boxy*y_height,...
                'FaceColor',cw,'FaceAlpha',0.2,...
                'LineWidth',1,'EdgeColor',cw);
            plot(zp(j)+(0+p)*dz,y0w,'Marker','+','Color',cw);
        end
    end
end