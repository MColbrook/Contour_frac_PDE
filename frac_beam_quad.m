function [zi,Vec]=frac_beam_quad(PDE,fhat,u0,u1,NUM,SING,OMEGA)
	opts = optimset('Display','none');
    if NUM.CONT==1 % algorithm 1 of paper
        h=lambertw(NUM.t1/NUM.t0*NUM.Nquad*pi*(pi-2*NUM.delta)/(NUM.beta*sin(pi/4-NUM.delta/2))*(1-sin((pi-2*NUM.delta)/4)))/NUM.Nquad;
        mu=NUM.beta/(NUM.t1*(1-sin(pi/4-NUM.delta/2)));
        alpha=pi/4-NUM.delta/2+h*mu*NUM.t1/(4*pi);
        ui=(-NUM.Nquad:NUM.Nquad)*h;
        zi=mu*(1+sin(1i*ui-alpha))+NUM.shift;
        zi_prime=mu*1i*cos(1i*ui-alpha);
    elseif NUM.CONT==2
        h=lambertw(NUM.t1/NUM.t0*NUM.Nquad*pi*(pi-2*NUM.delta)/(NUM.beta*sin(pi/4-NUM.delta/2))*(1-sin((pi-2*NUM.delta)/4)))/NUM.Nquad;
        mu=NUM.beta/(NUM.t1*(1-sin(pi/4-NUM.delta/2)));
        alpha=pi/4-NUM.delta/2+h*mu*NUM.t1/(4*pi);
        fun = @(X) max([-2*pi*(pi/2-X(3)-NUM.delta)/X(1),X(2)*NUM.t1-2*pi*X(3)/X(1),X(2)*NUM.t0*(1-sin(X(3))*cosh(NUM.Nquad*X(1))),...
            NUM.t1*(X(2)*(1-sin(X(3)))+NUM.shift)+log(NUM.tol/10)]);
        X=fminsearch(fun,[h,mu,alpha],opts);
        mu=X(2); h=X(1); alpha=X(3);
        ui=(-NUM.Nquad:NUM.Nquad)*h;
        zi=mu*(1+sin(1i*ui-alpha))+NUM.shift;
        zi_prime=mu*1i*cos(1i*ui-alpha);
    else % algorithm 2 of paper
        mu=NUM.beta/NUM.t1+1/(4*NUM.bb);
        h=(NUM.Nquad)^(-2/3)*(2*pi/(mu*NUM.t0)*(1-1/(2*sqrt(mu*NUM.bb))))^(1/3);
        fun = @(X) max([-2*pi/X(2)*(1-1/(2*sqrt(X(1)*NUM.bb))),-NUM.t1/(4*NUM.bb)-pi^2/(X(1)*NUM.t1*X(2)^2)+2*pi/X(2),-NUM.t0/(4*NUM.bb)-pi^2/(X(1)*NUM.t0*X(2)^2)+2*pi/X(2),...
            -NUM.t0/(4*NUM.bb)+X(1)*NUM.t0*(1-(X(2)*NUM.Nquad)^2),-NUM.t1/(4*NUM.bb)+X(1)*NUM.t1*(1-(X(2)*NUM.Nquad)^2),NUM.t1*(X(1)-1/(4*NUM.bb))+log(NUM.tol/10)]);
        X=fminsearch(fun,[mu,h],opts);
        mu=X(1); h=X(2);
        ui=(-NUM.Nquad:NUM.Nquad)*h;
        sigma=-1/(4*NUM.bb);
        zi=sigma+mu*(1i*ui+1).^2;
        zi_prime=2i*mu*(1i*ui+1);
    end
    
    III=find(abs(h*zi_prime.*exp(zi*NUM.t0))>10^(-16)); % get rid of negligible quadrature points
    zi=zi(III);
    zi_prime=zi_prime(III);
    NUM.Nquad=round((length(zi)-1)/2);
    
    % EXTRA CODE FOR SECTION 3
    if abs(imag(zi(abs(real(zi))==min(abs(real(zi))))))<OMEGA
        SING.CAUCHY=1;
    else
        SING.CAUCHY=0;
    end
    
    if NUM.Adaptive=="on" % adaptively select the discretisation size
        N=NUM.DiscMin;
        errIndx=1:length(zi);
        vi=cell(1,length(zi));
        LEN=zeros(length(zi),1)+N; % to keep track of lengths of expansions
        
        if NUM.Prog=="on"
            fprintf('Number of quadrature points = %d, Progress:',length(zi))
            pf = parfor_progress(length(zi));
            pfcleanup = onCleanup(@() delete(pf));
        end

        while(sum(errIndx)>0 && N<=NUM.DiscMax)
            [A,B,T_bc,S04,U,V]=elastic_beam_mats(PDE.a,PDE.b,PDE.rho_inv,N,PDE.BC);
            u0_c=chebcoeffs(u0,N,'kind',1);
            u1_c=chebcoeffs(u1,N,'kind',1);
            
            Z_temp=zi(errIndx); Z_prime_temp=zi_prime(errIndx); V_temp=vi(1,errIndx); temp=vi(1,errIndx);
            E_est=zeros(length(Z_temp),1);
            warning('off','all')
           
            if NUM.Parallel=="off"
                for j=1:length(Z_temp)
                    F_c=chebcoeffs(fhat(Z_temp(j))*PDE.rho_inv,N,'kind',1); % rescaled (divide by rho) forcing
                    if PDE.FRAC_TYPE=='RL'
                        if PDE.nu==1
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c);
                        else
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c);
                        end
                    elseif PDE.FRAC_TYPE=='C'
                        if PDE.nu>1
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c+Z_temp(j)^(PDE.nu-2)*u1_c);
                        else
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c);
                        end
                    end
                    RHS=[sparse(4,1);RHS(1:N-4)];

                    T=[T_bc;Z_temp(j)^2*S04(1:N-4,:)+A(1:N-4,:)+Z_temp(j)^PDE.nu*B(1:N-4,:)];

                    v=woodbury_solve(T,U,V,RHS);
                    V_temp{j}=v*h/(2*pi*1i)*Z_prime_temp(j);
                    E_est(j)=norm([temp{j}(:);zeros(N-length(temp{j}),1)]-V_temp{j}(:));
                end
            else
                parfor j=1:length(Z_temp)
                    F_c=chebcoeffs(fhat(Z_temp(j))*PDE.rho_inv,N,'kind',1); % rescaled (divide by rho) forcing
                    RHS=0;
                    if PDE.FRAC_TYPE=='RL'
                        if PDE.nu==1
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c);
                        else
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c);
                        end
                    elseif PDE.FRAC_TYPE=='C'
                        if PDE.nu>1
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c+Z_temp(j)^(PDE.nu-2)*u1_c);
                        else
                            RHS=S04*(F_c+Z_temp(j)*u0_c+u1_c)+B*(Z_temp(j)^(PDE.nu-1)*u0_c);
                        end
                    end
                    RHS=[sparse(4,1);RHS(1:N-4)];

                    T=[T_bc;Z_temp(j)^2*S04(1:N-4,:)+A(1:N-4,:)+Z_temp(j)^PDE.nu*B(1:N-4,:)];

                    v=woodbury_solve(T,U,V,RHS);
                    V_temp{j}=v*h/(2*pi*1i)*Z_prime_temp(j);
                    E_est(j)=norm([temp{j}(:);zeros(N-length(temp{j}),1)]-V_temp{j}(:));
                end
            end
            vi(errIndx)=V_temp;
            LEN(errIndx)=N;
            errIndx=find(E_est>max(10^(-14),NUM.tol/(max(1,sum(exp(max([NUM.t1*real(zi(:)),NUM.t0*real(zi(:))],[],2)))))));
            if NUM.Prog=="on"
                for j=1:(length(E_est)-length(errIndx))
                    parfor_progress(pf);
                end
            end
            N=round(2*N);
            warning('on','all')
        end
        
        Vec=zeros(max(LEN),length(zi));
        for j=1:length(zi)
            Vec(:,j)=[vi{j};zeros(max(LEN)-LEN(j),1)];
        end
    else
        N=NUM.N;
        Vec=zeros(N,length(zi));
    
        [A,B,T_bc,S04,U,V]=elastic_beam_mats(PDE.a,PDE.b,PDE.rho_inv,N,PDE.BC);
        u0_c=chebcoeffs(u0,N,'kind',1);
        u1_c=chebcoeffs(u1,N,'kind',1);
        
        if NUM.Prog=="on"
            fprintf('Number of quadrature points = %d, Progress:',length(zi))
            pf = parfor_progress(length(zi));
            pfcleanup = onCleanup(@() delete(pf));
        else
            pf=[];
        end

        if NUM.Parallel=="off"
            for j=1:length(zi)
                F_c=chebcoeffs(fhat(zi(j))*PDE.rho_inv,N,'kind',1); % rescaled (divide by rho) forcing
                if PDE.FRAC_TYPE=='RL'
                    if PDE.nu==1
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c);
                    else
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c);
                    end
                elseif PDE.FRAC_TYPE=='C'
                    if PDE.nu>1
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c+zi(j)^(PDE.nu-2)*u1_c);
                    else
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c);
                    end
                end
                RHS=[sparse(4,1);RHS(1:N-4)];

                T=[T_bc;zi(j)^2*S04(1:N-4,:)+A(1:N-4,:)+zi(j)^PDE.nu*B(1:N-4,:)];
                
                v=woodbury_solve(T,U,V,RHS);
                Vec(:,j)=v*h/(2*pi*1i)*zi_prime(j);
                if NUM.Prog=="on" % display progress
                    parfor_progress(pf);
                end
            end
        else
            parfor j=1:length(zi)
                F_c=chebcoeffs(fhat(zi(j))*PDE.rho_inv,N,'kind',1); % rescaled (divide by rho) forcing
                RHS=0;
                if PDE.FRAC_TYPE=='RL'
                    if PDE.nu==1
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c);
                    else
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c);
                    end
                elseif PDE.FRAC_TYPE=='C'
                    if PDE.nu>1
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c+zi(j)^(PDE.nu-2)*u1_c);
                    else
                        RHS=S04*(F_c+zi(j)*u0_c+u1_c)+B*(zi(j)^(PDE.nu-1)*u0_c);
                    end
                end
                RHS=[sparse(4,1);RHS(1:N-4)];

                T=[T_bc;zi(j)^2*S04(1:N-4,:)+A(1:N-4,:)+zi(j)^PDE.nu*B(1:N-4,:)];
                
                v=woodbury_solve(T,U,V,RHS);
                Vec(:,j)=v*h/(2*pi*1i)*zi_prime(j);
                if NUM.Prog=="on" % display progress
                    parfor_progress(pf);
                end
            end
        end
    end
    
    % now deal with poles of fhat
    if SING.CAUCHY==1
        if NUM.Adaptive=="on" % adaptively select the discretisation size
            N=NUM.DiscMin; errIndx=1:length(SING.POLES); vi_sing=cell(1,length(SING.POLES));
            LEN=zeros(length(SING.POLES),1)+N; % to keep track of lengths of expansions
            while(sum(errIndx)>0 && N<=NUM.DiscMax)
                warning('off','all')
                [A,B,T_bc,S04,U,V]=elastic_beam_mats(PDE.a,PDE.b,PDE.rho_inv,N,PDE.BC);
                Z_temp=SING.POLES(errIndx); V_temp=vi_sing(1,errIndx); temp=vi_sing(1,errIndx); E_est=zeros(length(Z_temp),1);
                if NUM.Parallel=="off"
                    for j=1:length(Z_temp)
                        RHS=S04*chebcoeffs(SING.RES{errIndx(j)}(Z_temp(j))*PDE.rho_inv,N,'kind',1);    RHS=[sparse(4,1); RHS(1:N-4)];
                        T=[T_bc;Z_temp(j)^2*S04(1:N-4,:)+A(1:N-4,:)+Z_temp(j)^PDE.nu*B(1:N-4,:)];
                        V_temp{j}=woodbury_solve(T,U,V,RHS);
                        E_est(j)=norm([temp{j}(:);zeros(N-length(temp{j}),1)]-V_temp{j}(:));
                    end
                else
                    parfor j=1:length(Z_temp)
                        RHS=S04*chebcoeffs(SING.RES{errIndx(j)}(Z_temp(j))*PDE.rho_inv,N,'kind',1);    RHS=[sparse(4,1); RHS(1:N-4)];
                        T=[T_bc;Z_temp(j)^2*S04(1:N-4,:)+A(1:N-4,:)+Z_temp(j)^PDE.nu*B(1:N-4,:)];
                        V_temp{j}=woodbury_solve(T,U,V,RHS);
                        E_est(j)=norm([temp{j}(:);zeros(N-length(temp{j}),1)]-V_temp{j}(:));
                    end
                end
                vi_sing(errIndx)=V_temp;
                LEN(errIndx)=N;
                errIndx=find(E_est>max(10^(-14),NUM.tol/(max([max(abs(exp(NUM.t0*SING.POLES))),max(abs(exp(NUM.t1*SING.POLES)))]))));
                N=round(2*N);
                warning('on','all')
            end
            Vec_sing=zeros(max(LEN),length(SING.POLES));
            for j=1:length(SING.POLES)
                Vec_sing(:,j)=[vi_sing{j};zeros(max(LEN)-LEN(j),1)];
            end
        else
            N=NUM.N;
            Vec_sing=zeros(N,length(SING.POLES));
            [A,B,T_bc,S04,U,V]=elastic_beam_mats(PDE.a,PDE.b,PDE.rho_inv,N,PDE.BC);
            if NUM.Parallel=="off"
                for j=1:length(SING.POLES)
                    RHS=S04*chebcoeffs(SING.RES{j}(SING.POLES(j))*PDE.rho_inv,N,'kind',1);    RHS=[sparse(4,1); RHS(1:N-4)];
                    T=[T_bc;SING.POLES(j)^2*S04(1:N-4,:)+A(1:N-4,:)+SING.POLES(j)^PDE.nu*B(1:N-4,:)];
                    Vec_sing(:,j)=woodbury_solve(T,U,V,RHS);
                end
            else
                parfor j=1:length(SING)
                    RHS=S04*chebcoeffs(SING.RES{j}(SING.POLES(j))*PDE.rho_inv,N,'kind',1);    RHS=[sparse(4,1); RHS(1:N-4)];
                    T=[T_bc;SING.POLES(j)^2*S04(1:N-4,:)+A(1:N-4,:)+SING.POLES(j)^PDE.nu*B(1:N-4,:)];
                    Vec_sing(:,j)=woodbury_solve(T,U,V,RHS);
                end
            end
        end
        
        zi=[zi,SING.POLES];
        L1=size(Vec,1); L2=size(Vec_sing,1);
        if L1>L2
            Vec=[Vec,[Vec_sing;zeros(L1-L2,length(SING.POLES))]];
        else
            Vec=[[Vec;zeros(L2-L1,size(Vec,2))],Vec_sing];
        end
    end
    
            
end

function [A,B,T_bc,S04,U,V]=elastic_beam_mats(a,b,rho_inv,N,BC)
    % computes the matrices corresponding to ultraspherical spectral method
    a_x=diff(a);    a_xx=diff(a_x);
    b_x=diff(b);    b_xx=diff(b_x);

    Mrho=ultraS.multmat(N,rho_inv,4);
    Ma=ultraS.multmat(N,a,4); Ma_x=ultraS.multmat(N,a_x,4); Ma_xx=ultraS.multmat(N,a_xx,4);
    Mb=ultraS.multmat(N,b,4); Mb_x=ultraS.multmat(N,b_x,4); Mb_xx=ultraS.multmat(N,b_xx,4);

    D2=ultraS.diffmat(N,2); D3=ultraS.diffmat(N,3); D4=ultraS.diffmat(N,4);

    S04=ultraS.convertmat(N,0,3); S24=ultraS.convertmat(N,2,3); S34=ultraS.convertmat(N,3,3);

    A=Mrho*(Ma*D4+2*Ma_x*S34*D3+Ma_xx*S24*D2);
    B=Mrho*(Mb*D4+2*Mb_x*S34*D3+Mb_xx*S24*D2);
    
    if BC==1 %1 is CC, 2 is SS, 3 is CS, 4 is SC (C=clamped, S=supported)
        V=[ones(1,N);  (-1).^(0:N-1);  (0:N-1).^2;  (-1).^(1:N).*(0:N-1).^2];
    elseif BC==2
        V=[ones(1,N);  (-1).^(0:N-1);  (((0:N-1).^4)-((0:N-1).^2))/3;  (-1).^(0:N-1).*(((0:N-1).^4)-((0:N-1).^2))/3];
    elseif BC==3
        V=[ones(1,N);  (-1).^(0:N-1);  (((0:N-1).^4)-((0:N-1).^2))/3;  (-1).^(1:N).*(0:N-1).^2];
    else
        V=[ones(1,N);  (-1).^(0:N-1);  (0:N-1).^2;  (-1).^(0:N-1).*(((0:N-1).^4)-((0:N-1).^2))/3];
    end
    T_bc=sparse(4,N);
    T_bc(1,1)=V(1,1); T_bc(2,2)=V(2,2); T_bc(3,3)=V(3,3); T_bc(4,4)=V(4,4);
    V(1,1)=0; V(2,2)=0; V(3,3)=0; V(4,4)=0; 
    U=[speye(4);sparse(N-4,4)];
end


function [x]=woodbury_solve(A,U,V,b)
    % solves (A+UV)x=b using Woodbury matrix identity
    M=V*(A\U);
    M=M+speye(size(M));
    x1=A\b;
    x=x1-A\(U*(M\(V*x1)));
end
