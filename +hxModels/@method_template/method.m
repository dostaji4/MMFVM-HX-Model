classdef method < hxModels.hxModelCommon
    
    properties (Access = private)        
    end
    
%% Implementation of Abstract methods    
    methods (Access = {?hxModel,?hxModels.hxModelCommon})
        % y = [Q, Two]
        % u = [m_dot, V_dot, Twi, Tai]
        [y, t] = step(obj,t, u, Ts);
    end
    
    methods (Access = protected)
        function initialize_hx(obj)             
           
        end
    end
    
    methods     
        
% Performance tips: 
% https://blogs.mathworks.com/loren/2012/03/26/considering-performance-in-object-oriented-matlab-code/

%% Method specific implementation (overloading)
        
%         function plot_outputs(obj,varargin)%            
%         end
        
%         function plot_states(obj,varargin)            
%         end
        
%         [Q, Two, t_y] = run(obj, m_dot, V_dot, Twi, Tai, Ts_u)

        % y = [Q, Two]
        % u = [m_dot, V_dot, Twi, Tai]        
%         [yss, xss] = ss(obj,u);         

    end

    methods(Static)
        
    end
end

