function Zeta = int_matrix_theta(n,K,M,L)
	Gamma = int_matrix(K,n,L);   
	Psi = selector_matrix(n,K,M);
    Zeta = [Gamma zeros((K-1)*(n+1),M); Psi eye(M)];
    function Psi = selector_matrix(n,K,M)
        Psi = zeros(M,(K-1)*(n+2));
        Psi(n+1,:)= ones(1,(K-1)*(n+2));
    end
end