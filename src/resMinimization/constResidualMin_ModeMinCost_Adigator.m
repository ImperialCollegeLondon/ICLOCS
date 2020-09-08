function [ const ] = constResidualMin_ModeMinCost_Adigator( X, U, P, T, data, mode)

X_s=X;
U_s=U;
P_s=P;

if isfield(data.dataNLP.data,'Xscale')
    X=scale_variables_back( X, data.dataNLP.data.Xscale_back, data.dataNLP.data.Xshift );
    U=scale_variables_back( U, data.dataNLP.data.Uscale_back, data.dataNLP.data.Ushift );
    if isfield(data.dataNLP.data,'Pscale')
        P=scale_variables_back( P, data.dataNLP.data.Pscale_back, data.dataNLP.data.Pshift );
    end
end


dataNLP=data.dataNLP;

t0=T(1);tf=T(end);
x0=X(1,:);
u0=U(1,:)';
xf=X(end,:)';
uf=U(end,:)';

%%

[~,~,f,g,avrc,b]=deal(data.dataNLP.functions_unscaled{:});

switch dataNLP.options.discretization
    
    case{'globalLGR','hpLGR'} % p/hp Transcription Method
        [~,~,~,~,ng,~,M,~,~,~,~,~,~,~,~,~,~]=deal(dataNLP.sizes{1:17});
        
        if mode==1
            X_mesh=X(1:end-1,:);
            U_mesh=U(1:end-1,:);
            T_mesh=T(1:end-1,:);
            X_mesh_Np1=X;

            g_vect=reshape(g(X_mesh,U_mesh,P,T_mesh,data.dataNLP.data)',M*ng,1);
            if isfield(data.dataNLP.data,'Xscale')
                cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
            else
                cr=avrc(X_mesh_Np1,U_mesh,P,T,data.dataNLP)';
            end
            const=[  g_vect(data.dataNLP.gAllidx);
                     cr;
                     b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{})];

        else
            const=[];
        end
                
    otherwise
        [~,~,n,~,ng,~,M,~,~,~,~,~,~]=deal(dataNLP.sizes{1:13});
        
        if mode==1
            g_vect=reshape(g(X,U,P,T,data.dataNLP.data)',M*ng,1);
            if isfield(data.dataNLP.data,'Xscale')
                cr=avrc(X_s,U_s,P_s,T,data.dataNLP)';
            else
                cr=avrc(X,U,P,T,data.dataNLP)';
            end
            const=[ g_vect(data.dataNLP.gAllidx);
                    cr;
                    b(x0,xf,u0,uf,P,t0,tf,data.dataNLP.data,{})];
        else
            const=[];
        end
end
      
end

