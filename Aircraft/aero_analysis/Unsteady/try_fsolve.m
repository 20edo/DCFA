% Following iterations
for i=2:length(v)
    X = zeros(size(X_old));
    e = zeros(size(e_old));
    for k=1:size(X_old,2)
        tic
        res = 1;
        it=0;
        while res > 1e-3 && it<10
            if it ==0
                % Build matrices to solve the non-linear system
                kk = l*imag(e_old(k))/v(i);
                wing = m_add_unsteady_loads(wing,[1,0,0]',kk);
                Ham = wing.Ham;
                Ham_dk = wing.Ham_dk;
                Ham = V'*Ham*V;
                Ham_dk = V'*Ham_dk*V;
                A(1:size(M,1),1:size(M,1))=e_old(k)^2*M+e_old(k)*Cs+K-q(i)*(Ham);
                A(1:size(M,1),end)=(2*e_old(k)*M+Cs-q(i)*(-1i*Ham_dk*l/v(i)))*X_old(:,k);
                A(end,1:size(M,1))=2*X_old(:,k)';
                A(end,end)=0;
                b(1:size(M,1),1)=-A(1:size(M,1),1:size(M,1))*X_old(:,k);
                b(end)=1-X_old(:,k)'*X_old(:,k);
            end
            %         funz=@(z) A*z-b;
            %         z0=[X_old(:,k);e_old(k)];
            %         [z,~,exitflag]=fsolve(funz,z0);
            z=A\b;
            X(:,k)=X_old(:,k)+z(1:end-1);
            e(k)=e_old(k)+z(end);
            phrase = ['Eig number ',num2str(k),' out of ',num2str(length(e_old)),'; Velocity ',num2str(i),' out of ',num2str(length(v))];
            % Build matrices to calculate residual
            kk = l*imag(e_old(k))/v(i);
            wing = m_add_unsteady_loads(wing,[1,0,0]',kk);
            Ham = wing.Ham;
            Ham_dk = wing.Ham_dk;
            Ham = V'*Ham*V;
            Ham_dk = V'*Ham_dk*V;
            A(1:size(M,1),1:size(M,1))=e_old(k)^2*M+e_old(k)*Cs+K-q(i)*(Ham);
            A(1:size(M,1),end)=(2*e_old(k)*M+Cs-q(i)*(-1i*Ham_dk*l/v(i)))*X_old(:,k);
            A(end,1:size(M,1))=2*X_old(:,k)';
            A(end,end)=0;
            b(1:size(M,1),1)=-A(1:size(M,1),1:size(M,1))*X_old(:,k);
            b(end)=1-X_old(:,k)'*X_old(:,k);
            % Calculate residual and update it
            res = norm((e_old(k)^2*M+e_old(k)*Cs+K-q(i)*(Ham))*X_old(:,k));
            it= it+1;
        end
        disp(phrase)
        toc
        
    end
    X_old = X;
    e_old = e;
    eig_(i,:) = e_old;
end