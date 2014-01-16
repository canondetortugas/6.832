function [K, S] = lqr_stabilizer(A, B, Q, R, dt)

N = length(A);

S = cell(1,N);
K = cell(1,N);

[Kf, Sf] = lqr(A{end}, B{end}, Q, R);
S{end} = Sf;
K{end} = Kf;

for idx = (N-1):-1:1
    
    Sp = S{idx+1};
    Ap = A{idx+1};
    Bp = B{idx+1};

    S{idx} =S{idx+1} + (Q - Sp*Bp*inv(R)*Bp'*Sp + Sp*Ap + Ap'*Sp).*dt;
    K{idx} = inv(R)*B{idx}'*S{idx};

end


end