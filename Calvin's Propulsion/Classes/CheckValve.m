classdef CheckValve < matlab.mixin.Copyable
    properties
        name
        type
        Cv
        p_crack
        r_pop
        phase
    end
    methods
        function obj = CheckValve()
            obj.phase = 'gas';
            obj.type = 'valve';
        end
        function q = CheckValveGasFlowRate(obj,p_1,p_2,T_1,gas)
            dp = p_1 - p_2;
            if (dp < obj.p_crack)
                Cv_eff = 0;
            elseif (dp - obj.p_crack < obj.p_crack*obj.r_pop)
                Cv_eff = obj.Cv*sqrt((dp - obj.p_crack)/(obj.p_crack*obj.r_pop));
            else
                Cv_eff = obj.Cv;
            end
            q = ValveGasFlow(Cv_eff,p_1,p_2,T_1,gas.gam,gas.G);
        end
        function mdot = GasMassFlowRate(obj,p_1,p_2,T_1,gas)
            q = CheckValveGasFlowRate(obj,p_1,p_2,T_1,gas);
            mdot = q*14.7*144/(gas.R*528)/60;
        end
    end
end