classdef Gao_Polyhedron
    % a Polyhedron defined by dengwei gao
    properties
        V
        Dim
    end
    
    methods
        function obj = Gao_Polyhedron(V)
            obj.V = V;
            obj.Dim = size(V,2);
        end
        
        function obj = mtimes(T,obj)
            obj.V = obj.V*T';
        end
        
        function obj = plus(obj,S)
            Vp = obj.V; Vs = S.V; 
            
			if size(Vp, 1)==0
				obj.V = Vs;
			elseif size(Vs, 1)==0
				obj.V = Vp;
			else
				obj.V = zeros(size(Vp,1)*size(Vs,1), obj.Dim);
				for i = 1:size(Vs, 1)
					obj.V((i-1)*size(Vp, 1)+1:i*size(Vp, 1), :) = bsxfun(@plus, Vp, Vs(i, :));
				end
			end
            ind = convhull(obj.V);  % only need the covex hull
            obj.V = obj.V(ind,:);
        end
    end
end

