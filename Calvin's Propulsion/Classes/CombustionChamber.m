classdef CombustionChamber < matlab.mixin.Copyable
    properties
        name
        type
        cT
        p
        cstar
        OF
        mdot
        FT
        A_t
        data
    end
    methods
        function obj = CombustionChamber()
            obj.type = 'volume';
            obj.p = 0;
            obj.cstar = 0;
            obj.FT = 0;
            obj.data = struct('p',[],'OF',[],'cstar',[],'FT',[]);
        end
        function obj = Thrust(obj)
            if abs(obj.p*obj.A_t*obj.cT - obj.mdot*obj.cstar*obj.cT) > 1
                error('Thrust error')
            else
                obj.FT = obj.p*obj.A_t*obj.cT;
            end
        end
        function obj = ZeroDataArrays(obj,n_sim)
            obj.data.p = zeros(1,n_sim);
            obj.data.OF = zeros(1,n_sim);
            obj.data.cstar = zeros(1,n_sim);
            obj.data.FT = zeros(1,n_sim);
        end
        function obj = RecordStates(obj,i)
            obj.data.p(i) = obj.p;
            obj.data.OF(i) = obj.OF;
            obj.data.cstar(i) = obj.cstar;
            obj.data.FT(i) = obj.FT;
        end
    end
end