classdef Regulator < matlab.mixin.Copyable
    properties
        name
        type
        Cv
        p_set
        p_crack
        r_droop
        p_droop_cali
        SPE
        phase
    end
    methods
        function obj = Regulator()
            obj.phase = 'gas';
            obj.type = 'valve';
        end
        
        function mdot = GasMassFlowRate(obj,loc,map,Xvals,A,dt,gas)
            
            [h_reg,w_reg] = matsplit(loc);
            
            p_1 = TargetVolumeParameter(map,[h_reg,w_reg],'in','p');
            T_1 = TargetVolumeParameter(map,[h_reg,w_reg],'in','T');
            p_2 = TargetVolumeParameter(map,[h_reg,w_reg],'out','p');
            
            q_max = ValveGasFlow(obj.Cv,p_1,p_2,T_1,gas.G,gas.gam);
                
            n = 10;
            nrange = linspace(1,n,n);
            mdotsv = cell(1,10);

            flowrange = linspace(0,q_max,n);
            mdotrange = flowrange*14.7*144/(gas.R*528)/60;

            for v = 1:1:n
                newmap = CopyMap(map);
                newmap = DealGasMassFlowRates(newmap,[h_reg,w_reg],mdotrange(v));
                newvolumes = VectorizeMapVolumes(newmap);
                
                for j = 1:1:length(newvolumes)
                    newvolumes{j}.summdots;
                    mdotsv{v} = vertcat(mdotsv{v},newvolumes{j}.vmdot);
                end
            end
            
            pseudomap = CopyMap(map);
            oldXvals = Xvals;
            pseudop_1 = zeros(1,n);
            pseudoT_1 = zeros(1,n);
            pseudop_2 = zeros(1,n);
            p_effrange = zeros(1,n);
            for v = 1:1:n
                pseudovolumes = VectorizeMapVolumes(pseudomap);
                pseudomdots = mdotsv{v};
                dX = double(A\pseudomdots);
                newXvals = oldXvals + dX*dt;

                DealXvals(pseudovolumes,newXvals);
                
                pseudop_1(v) = TargetVolumeParameter(pseudomap,[h_reg,w_reg],'in','p');
                pseudoT_1(v) = TargetVolumeParameter(pseudomap,[h_reg,w_reg],'in','T');
                pseudop_2(v) = TargetVolumeParameter(pseudomap,[h_reg,w_reg],'out','p');
                
            end
            for v = 1:1:n
                p_effrange(v) = DroopPressure(obj,flowrange(v),p_1,p_2,T_1,gas);
            end
            q_targ = 0;
            if (max(p_effrange) < min(pseudop_2))
            else
                for v = 1:1:n-1
                    if (p_effrange(v) > pseudop_2(v)) && (p_effrange(v+1) < pseudop_2(v+1))
                        n_targ = v + (pseudop_2(v) - p_effrange(v))/(p_effrange(v+1) - p_effrange(v) - pseudop_2(v+1) + pseudop_2(v));
                        q_targ = lininterp1(nrange,flowrange,n_targ);
                        break
                    end
                end
                if q_targ == 0 && p_1 > p_2
                    q_targ = ValveGasFlow(obj.Cv,p_1,p_2,T_1,gas.G,gas.gam);
                end
            end
            mdot = q_targ*14.7*144/(gas.R*528)/60;     
        end
        function p_eff = DroopPressure(obj,q,p_1,p_2,T_1,gas)
            q_max = ValveGasFlow(obj.Cv,p_1,p_2,T_1,gas.G,gas.gam);
            if q <= q_max
                p_eff = obj.r_droop*obj.p_droop_cali/p_1*q*sqrt(gas.G) + obj.p_set - obj.p_crack - obj.SPE*(p_1 - obj.p_set); %(q/q_max)
            end 
        end
    end
end