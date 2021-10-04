classdef FVM < hxModels.hxModelCommon
    properties 
       flux_limiter = 'HCUS';  % one of limiterList        
       integration_method = 'RK2';
       simulation_mode = 'discrete'; %'continuous';
    end
    
    properties (Constant)
        model_type = 'FVM';
    end
%         
    properties (Access = {?hxModel, ?hxModelCommon}, Constant)
        flux_limiterList ={'Van_Leer', 'superbee', 'smart', 'mc', 'minmod',...
              'Koren', 'CHARM', 'Osher', 'UMIST','ospre','HCUS','HQUICK',...
              'van_Albada_1','van_Albada_2','mg','Sweby','fo'};  
        integration_methodList = {'RK1','RK2','RK2e','RK4'}        
    end
    properties (Access = private)        
        iLimiter = 0; %Limiter index
        iIntMethod = 0; %Integration mode index                      
    end
    
    %% Method specific implementation
    methods
        
    end
    
    %% Prototypes
    methods
        % y = [Q, Two]
        % u = [m_dot, V_dot, Twi, Tai]
        [Q, Two, t_y] = run(obj, m_dot, V_dot, Twi, Tai, Ts_u);        
    end
    
    methods (Access = {?hxModel,?hxModels.hxModelCommon})
        % prototype
        [y, t] = step(obj,t, u, Ts)       
        dx = dxdt_rhs(obj,t,x,u,Ts_u,z,iLim);       
    end
    
    methods (Access = protected)

        function initialize_hx(obj)        
            
            % --- Method specific assertions
            props = {'simulation_mode','flux_limiter','integration_method'};
            hasMemberList = [1,1,1];
            for i = 1:numel(props)
                if hasMemberList(i) 
                    % Check non-emptiness
                    assert(~isempty(obj.(props{i})),...
                        'Property "%s" not specified. Set the property with one of the following options: ''%s''.',...
                        props{i}, strjoin(obj.([props{i} 'List']),''', '''));
                    % Check validity                
                    assert(ismember(obj.(props{i}),obj.([props{i} 'List'])),...
                        'Unknown setting "%s" for property "%s". Use one the the following options: ''%s''.',...
                        obj.(props{i}),props{i},strjoin(obj.([props{i} 'List']),''', '''));
                else
                    % Check non-emptiness
                    assert(~isempty(obj.(props{i})), 'Property "%s" not set. Please specify it first.',...
                        props{i});
                end
            end
            obj.iLimiter = find(ismember(obj.flux_limiterList,obj.flux_limiter)); % for performance reasons are limiters called by their index            
            obj.iIntMethod = find(ismember(obj.integration_methodList,obj.integration_method));
            
            % --- Initialize precomputed constants
            obj.Cw = obj.mw * obj.cp_w;  %[J/K] heat capacity of water inside HX
            obj.mw_e = obj.mw / obj.Nw;
%             if isempty(obj.Ts_y_max)
%                 obj.Ts_y_max = obj.mw_e / obj.m_dot_min; % Sampling interval used for energy flow model [s]
%             end
            
            % --- Initial water and body temperatures
            tmp = {obj.Tw0,obj.Tb0};
            tmpText = {'Tw0','Tb0'};
            nData = [obj.Nw, obj.Nb];
            for i = 1:numel(tmp)
                assert(~isempty(tmp{i}),'Initialize water and body temperatures first. Allowed values are: scalar for constant profile or a vector of dimension Nw=%d (Nb=%d).',obj.Nw,obj.Nb);
                if numel(tmp{i}) == 1
                    % Replicate nData times
                    tmp{i} = tmp{i} * ones(nData(i),1);
                elseif numel(tmp{i}) == nData(i)
                    % Ensure column vector
                    tmp{i} = tmp{i}(:);
                else
                    error('Incorrect dimension of "%s". Specify either a scalar value or a vector of length %d.',tmpText{i},nData(i));
                end
            end
            [obj.Tw0,obj.Tb0] = tmp{:};
            
            % --- State space model
            obj.x = [obj.Tw0; obj.Tb0]; %Initial state            
            obj.nx = obj.Nw + obj.Nb; % Default for FVM
            obj.ny = 2;
            
            obj.idx_Tw = 1:obj.Nw;
            obj.idx_Tb = obj.Nw+1:obj.Nw+obj.Nb;
            
            % --- Define grid (grid without boundary points)
            [obj.z,obj.dz] = hxModels.hxModelCommon.equidistant_grid(obj.Nw);            
        end
            
    end
    
    methods(Static)       
    end
end

