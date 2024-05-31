function T = int_transform(der,K,M,L,knot_switch)
	if ~exist('knot_switch','var')
        knot_switch = 1;
	end
	m = M-der;
	T = int_matrix_theta(0,K,M,L);
	for i = 1:m-1
        T = T*int_matrix_theta(i,K,M,L);
    end
	if knot_switch==1
        T = T*knot__matrix(m,K,M);
    end
    
    function P = knot__matrix(n,K,M)
        P = zeros((K-1)*(n+1)+M,K);
        k=1;
        for ii = 0:(K-1)*(n+1)-1
            for j = 0:K-1
                if ii==0 && j==0
                    P(ii+1,j+1)=1;
                elseif ii==k*(n+1)-1 && j==k && k<=K
                    P(ii+1,j+1)=1;
                    k=k+1;
                end
            end
        end
    end    
end