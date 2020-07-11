classdef PropellantTank  < matlab.mixin.Copyable
    properties
        name
        type
        p
        V
        T
        rho
        m
        rho_prop
        m_prop
        dE
        A
        Xvals
        vmdotin
        vmdotliq
        vmdot
        data
    end
    methods
        function obj = PropellantTank()
            obj.vmdotin = zeros(6,1);
            obj.vmdotliq = zeros(6,1);
            obj.type = 'volume';
            obj.data = struct('p',[],'V',[],'T',[],'rho',[],'m',[],'m_prop',[]);
        end
        function obj = ZeroDataArrays(obj,n_sim)
            obj.data.p = zeros(1,n_sim);
            obj.data.V = zeros(1,n_sim);
            obj.data.T = zeros(1,n_sim);
            obj.data.rho = zeros(1,n_sim);
            obj.data.m = zeros(1,n_sim);
            obj.data.m_prop = zeros(1,n_sim);
        end
        function obj = RecordStates(obj,i)
            obj.data.p(i) = obj.p;
            obj.data.V(i) = obj.V;
            obj.data.T(i) = obj.T;
            obj.data.rho(i) = obj.rho;
            obj.data.m_prop(i) = obj.m_prop;
        end
        function obj = DefineGoverningEquations(obj,upstreamobj,gas)
            ACoM1 = [0,0,-obj.V,1,-obj.rho,0];
            %CoM1 = (msym) - obj.V*(rhosym) - obj.rho*(Vsym) == 0;
            ACoM2 = [0,0,0,1,0,0];
            %CoM2 = (msym) == mdotinsym;
            AVC = [0,0,0,0,obj.rho_prop,0];
            %VC = obj.rho_prop*(Vsym) == mdot_liqsym;
            AEoS = [144,-gas.R*obj.rho,-gas.R*obj.T,0,0,0];
            %EoS = (psym)*144 - gas.R*obj.rho*(Tsym) - gas.R*obj.T*(rhosym) == 0;
            ACoE = [0,obj.m/(gas.gam*obj.T),0,obj.T/(gas.gam*upstreamobj.T),(gas.gam - 1)*obj.p*144/(gas.gam*gas.R*upstreamobj.T),0];
            %CoE = obj.m*(Tsym)/(gas.gam*obj.T) + obj.T*(msym)/(gas.gam*upstreamobj.T) + (gas.gam - 1)*obj.p*(Vsym)/(gas.gam*gas.R*upstreamobj.T) == mdotinsym;
            ACoMliq = [0,0,0,0,0,1];
            %CoMliq = (mliqsym) == -mdot_liqsym;
            %obj.Amdots = [0;mdotinsym;mdot_liqsym;0;mdotinsym;-mdot_liqsym];
            obj.A = [ACoM1; ACoM2; AVC; AEoS; ACoE; ACoMliq];
            
            if obj.m_prop <= 0
                notify(obj,'OutOfPropellant');
            end
        end
        function obj = DefineXvals(obj)
            obj.Xvals = [obj.p obj.T obj.rho obj.m obj.V obj.m_prop];
        end
        function obj = Resetmdots(obj)
            obj.vmdotin = zeros(6,1);
            obj.vmdotliq = zeros(6,1);
            obj.vmdot = zeros(6,1);
        end
        function obj = summdots(obj)
            obj.vmdot = obj.vmdotin + obj.vmdotliq;
        end
    end
    events
        OutOfPropellant
    end
end
    