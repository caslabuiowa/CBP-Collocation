function I = int_matrix(K,m,L)
        I = zeros((K-1)*(m+1),1);
        I(:,2:(K-1)*(m+1)+1) = (L/((K-1)*(m+1))).*triu(ones((K-1)*(m+1)));
        for i = 1:K-2
            I_front = I(:,1:(i)*(m+1)+i);
            I_back = I(:,(i)*(m+1)+i:end);
            I = [I_front, I_back];
        end
end