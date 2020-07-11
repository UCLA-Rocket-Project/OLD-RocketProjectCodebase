classdef PressurantTank  < matlab.mixin.Copyable
    properties
        name
        type
        p
        V
        T
        rho
        m
        A
        Xvals
        vmdotout
        vmdot
        data
    end
    methods
        function obj = PressurantTank()
            obj.vmdotout = zeros(4,1);
            obj.type = 'volume';
            obj.data = struct('p',[],'V',[],'T',[],'rho',[],'m',[]);
        end
        function obj = ZeroDataArrays(obj,n_sim)
            obj.data.p = zeros(1,n_sim);
            obj.data.V = zeros(1,n_sim);
            obj.data.T = zeros(1,n_sim);
            obj.data.rho = zeros(1,n_sim);
            obj.data.m = zeros(1,n_sim);
        end
        function obj = RecordStates(obj,i)
            obj.data.p(i) = obj.p;
            obj.data.V(i) = obj.V;
            obj.data.T(i) = obj.T;
            obj.data.rho(i) = obj.rho;
            obj.data.m(i) = obj.m;
        end
        function obj = DefineGoverningEquations(obj,gas)
            ACoM1 = [0,0,-obj.V,1];
            %CoM1 = (msym) - obj.V*(rhosym) == 0;
            ACoM2 = [0,0,0,1];
            %CoM2 = (msym) == -mdotouttotsym;
            AEoS = [144,-gas.R*obj.rho,-gas.R*obj.T,0];
            %EoS = (psym)*144 - gas.R*obj.rho*(Tsym) - gas.R*obj.T*(rhosym) == 0;
            ACoE = [0,obj.m/(gas.gam*obj.T),0,1/gas.gam];
            %CoE = obj.m*(Tsym)/(gas.gam*obj.T) + (msym)/gas.gam == -mdotouttotsym;
            obj.A = [ACoM1; ACoM2; AEoS; ACoE];
        end
        function obj = DefineXvals(obj)
            obj.Xvals = [obj.p obj.T obj.rho obj.m];
        end
        function obj = Resetmdots(obj)
            obj.vmdotout = zeros(4,1);
            obj.vmdot = zeros(4,1);
        end
        function obj = summdots(obj)
            obj.vmdot = obj.vmdotout;
        end
%         function obj = AssignXvals(obj,newXvals,Xvars)
%             for i = 1:1:length(newXvals)
%                 if Xvars(i) == obj.Xvars(1)
%                     obj.p = newXvals(i);
%                     obj.T = newXvals(i+1);
%                     obj.rho = newXvals(i+2);
%                     obj.m = newXvals(i+3);
%                     return;
%                 end
%             end
%         end
    end
end