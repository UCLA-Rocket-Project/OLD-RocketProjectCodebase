classdef InjectorManifold < matlab.mixin.Copyable
    properties
        name
        type
        CdA
        Cv
        phase
        side
    end
    methods
        function obj = InjectorManifold()
            obj.phase = 'liquid';
            obj.type = 'valve';
        end
        function obj = AssignType(side)
            if strcmp(side,'oxidizer') || strcmp(side,'fuel')
                obj.type = side;
            else
                error('injector manifold side not valid. must be fuel or oxidizer')
            end
        end
    end
end