function solution=multipleShooting(required,z,data)

% Generate the cost,constraint and gradient information needed for the
% multiple shooting method.
% Note: The multiple shooting segments are evaluated sequentially
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk


%--------------------------------------------------------------------------
% Initialization
%--------------------------------------------------------------------------

global sol;


if ~(isfield(sol,'z')) || length(z)~=length(sol.z) 
   sol.z=zeros(length(z),1);
end

% Assign variables from data structure

% Get dimensions of problem

[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});
nm=n+m;

% Extract p and tf
if nt; tf=z(1); else tf=data.tf; 
end
p=z(nt+1:nt+np);
vdat=data.data;
x0=z(nt+np+1:nt+np+n);


% Assign some data to a new structure that is passed to the ODE functions

params=data;
params.parameters=p;


% Options for Sundials
%data.options.cvodesf.SensRhsFn=@sensFn;


% Get handles for  functions needed to be evaluated in this scope
[L,E,f,g,avrc,b]=deal(data.functions{:});


% It is not possible to have tf as variable => nt=0; I need to introduce a errorcheck 


% Should it be t0???
%tk=0;

% It should be for the time varying case
t0=data.t0;
t=t0;
T=[0;cumsum(data.tau)];



%--------------------------------------------------------------------------
% Solve the IVP
%--------------------------------------------------------------------------


% If  NLP Solver supplies new variable solve the IVP
if (max(abs(z-sol.z))>eps)

%     % Generate random vector for jacobian structure determination if needed
%     if strcmp(required,'jacStruct')
%         z=rand(length(z),1);
%     end


    stageCost=0;            % Initialize stage cost
    c=(x0-data.x0t)*data.cx0;          % Consistency of initial conditions
    gz=[];                  % Path constraints
    
  
    jacgz=[];
    jaccz=[zeros(n,np), eye(n)*data.cx0, zeros(n,N*(n+m))];      %Jacobian of the conditions on the state dynamic
    
    grad=[];

    % Initial conditions for sensitivities
    ic_s=zeros(n+1,nm);ic_s(1:n,1:n)=eye(n);
    if ng % Add user-defined constraints if needed 
        gz=[gz;(g(z(1:n).',z(n+1:nm).',p',t0,vdat)).'];
        if data.FD.index.dg.flag
          [dfA,dgzA]=jacConst(f,g,z(1:n)',z(n+1:n+m)',p',t0,b,x0,z(1:n),z(n+1:n+m),z(n+1:n+m),p,tf,t0,data);
          dgz=reshape(cell2mat([dgzA.dx, dgzA.du]),ng,n+m);
          dgp=reshape(cell2mat([dgzA.dp]),ng,np);           
        else  
        params.input=z(n+1:n+m);   
        [dgp,dgz]=dgdz(z(1:n).',z(n+1:n+m).',p',t0,params);
        end
        jacgz=[jacgz;dgp, dgz, zeros(ng,N*(n+m)-m)];
    end

 
     
    for i=0:M-3 % Compute costs,constraints and sensitivities M-1 times
        k=i*nm+1;  
%       if strcmp(required,'jacStruct')
%          jac=[jac;zeros(n,k-1) ones(n,nm) -eye(n) zeros(n,(nm)*(M-i-2))];
%       else   
        params.input=z(k+n:k+nm-1);

        % Define the reference for the cost function in the i-th interval
        if ~isempty(data.references.xr)
          params.references.xr=data.references.xr(i+1,:);
        end
        if ~isempty(data.references.ur)
         params.references.ur=data.references.ur(i+1,:);
        end
        
        data.options.cvodes.UserData=params;
        
        % Allocate memory
        %CVodeMalloc(@odeFn,t,[z(k:k+n-1);0],data.options.cvodes,params);
        CVodeInit(@odeFn,data.options.cvodes.LMM,data.options.cvodes.NonlinearSolver,t,[z(k:k+n-1);0],data.options.cvodes);
        %CVodeSensMalloc(nm,'Staggered',ic_s,data.options.cvodesf);
        CVodeSensInit(nm,@sensFn,ic_s,data.options.cvodesf);
        
        % Solve the IVP, integrate stage cost and determine forward sensitivity
        % 
        [status,t,xf,xS] = CVode((tf-t0)*T(i+2)+t0,'Normal');

        % Free memory
        CVodeFree;

        % Store the solution in the required matrices. Update cost.
        stageCost=stageCost+xf(end);
        c=[c;z(k+nm:k+nm+n-1)-xf(1:n)];
        jaccz=[jaccz;zeros(n,k-1) -xS(1:n,:) eye(n) zeros(n,(nm)*(M-i-2))];
        grad=[grad xS(n+1,:)];

%       end 
    % Add user-defined path constraints if needed
        if ng
            gz=[gz;(g(z(k+nm:k+nm+n-1).',z(k+n:k+nm-1).',p',t,vdat)).'];
            if data.FD.index.dg.flag
             [dfA,dgzA]=jacConst(f,g,z(k+nm:k+nm+n-1)',z(k+n:k+nm-1)',p',t,b,x0,z(1:n),z(n+1:n+m),z(n+1:n+m),p,tf,t0,data);
             dgz=reshape(cell2mat([dgzA.dx, dgzA.du]),ng,n+m);
             dgp=reshape(cell2mat([dgzA.dp]),ng,np);           
            else    
              [dgp,dgz]=dgdz(z(k+nm:k+nm+n-1).',z(k+n:k+nm-1).',p',t,params);
            end
             jacgz=[jacgz;dgp, zeros(ng,k-1+nm),dgz,zeros(ng,nm*(M-i-2)-m)];
        end

    end
       

    k=(M-2)*nm+1;  
    params.input=z(k+n:k+nm-1);

    % Define the reference for the cost function in the i-th interval
      if ~isempty(data.references.xr)
        params.references.xr=data.references.xr(M-1,:);
      end
      if ~isempty(data.references.ur)
        params.references.ur=data.references.ur(M-1,:);
      end
       
        data.options.cvodes.UserData=params;
        
        % Allocate memory
        %CVodeMalloc(@odeFn,t,[z(k:k+n-1);0],data.options.cvodes,params);
        CVodeInit(@odeFn,data.options.cvodes.LMM,data.options.cvodes.NonlinearSolver,t,[z(k:k+n-1);0],data.options.cvodes);
        %CVodeSensMalloc(nm,'Staggered',ic_s,data.options.cvodesf);
        CVodeSensInit(nm,@sensFn,ic_s,data.options.cvodesf);
        
        

        % Solve the IVP, integrate stage cost and determine forward sensitivity
        % It is uniform (I need to generalize with tau)
        [status,t,xf,xS] = CVode((tf-t0)*T(M)+t0,'Normal');

        % Free memory
        CVodeFree;

        % Store the solution in the required matrices. Update cost.
        stageCost=stageCost+xf(end);
        c=[c;z(k+nm:k+nm+n-1)-xf(1:n)];
        jaccz=[jaccz;zeros(n,k-1) -xS(1:n,:) eye(n)];
        grad=[grad xS(n+1,:)];
    
    
    
    
  
    c=[c;gz];
    
    
  % Add user-defined terminal constraints if needed
    jacbz=[];
   
    if nb
        c=[c;b(z(1:n),z(end-n+1:end),z(n+1:nm),z(k+n:k+nm-1),p,0,tf,vdat)];
        jacbz=dbdz(z(1:n),z(end-n+1:end),z(n+1:nm),z(k+n:k+nm-1),p,tf,params);
    end


    jac=[jaccz;jacgz;jacbz];
     
    % Update the results and store in global structure zo

    sol.z=z;
    
%     stageCost
%     E(z(1:n),z(end-n+1:end),z(n+1:nm),z(k+n:k+nm-1),p,tf,vdat)
    sol.cost=stageCost+E(z(1:n),z(end-n+1:end),z(n+1:nm),z(k+n:k+nm-1),p,0,tf,vdat); 
    sol.const=c;
    sol.jacConst=sparse(jac);
   
    sol.gradCost=[grad zeros(1,n)]+dEdz(z,params);
    
end

%
%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------

switch required

    case{'cost'}
        
       solution=sol.cost;

    case{'gradCost'}
        
       solution=sol.gradCost;

    case{'const'}
        solution=sol.const;

    case{'jacConst'}
        
       solution=sol.jacConst;
   
end



%--------------------------------------------------------------------------
% Function definitons
% -------------------------------------------------------------------------

% Augmented ODE function with z:= [x;y]
% x'= f(x,u) and y'= L(x,u)

function [dz, flag, new_data] = odeFn(t, z, data)

u=data.input;
p=data.parameters;
f=data.functions{3};L=data.functions{1};
xr=data.references.xr;ur=data.references.ur;

vdat=data.data;


dz=[f(z(1:end-1)',u(:)',p',t,vdat)';L(z(1:end-1)',xr,u(:)',ur,p',t,vdat)]; 

flag = 0;new_data = [];

%------------------------------------------------------------------
% Forward sensitivity equations

function [dzS, flag, new_data] = sensFn(t,z,dz,zS,data)

u=data.input';
x=z(1:end-1)';
ur=data.references.ur;   
xr=data.references.xr;



ei=data.options.perturbation.J;
p=data.parameters';
[nt,np,n,m]=deal(data.sizes{1:4});
[L,E,f,g,avrc,b]=deal(data.functions{:});
vdat=data.data;

%  Compute dfdx: Check if analytic expression is available. If not
%  approximate with data.method


if data.FD.index.df.flag
    [df]=jacConst(f,g,x,u,p,t,b,x',x',data.input,data.input,p',0,0,data);
    dfdx_s=cell2mat(df.dx.').';
    dfdu_s=cell2mat(df.du.').';
else
     dfdx_s=[];
     for i=1:n
       e=zeros(1,n);e(i)=ei;
       dfdx_s=[dfdx_s (f(x+e,u,p,t,vdat)-f(x-e,u,p,t,vdat))'/(2*ei)]; % fix? ECK 15/9/08
     end
     dfdu_s=[];
     for i=1:m
        e=zeros(1,m);e(i)=ei;
        dfdu_s=[dfdu_s (f(x,u+e,p,t,vdat)-f(x,u-e,p,t,vdat))'/(2*ei)]; % fix? ECK 15/9/08
      end 
end


%  Compute dLdx: Check if analytic expression is available. If not
%  approximate with data.method

if data.FD.index.dL.flag
   [dL]=gradCost(L,x,xr,u,ur,p,t,E,x',x',data.input,data.input,p',t,data);
   dLdx_s=dL.dx;
   dLdu_s=dL.du;
else

            dLdx_s=[];
            for i=1:n
                e=zeros(1,n);e(i)=ei;
                dLdx_s=[dLdx_s (L(x+e,xr,u,ur,p,t,vdat)-L(x-e,xr,u,ur,p,t,vdat))/(2*ei)]; 
            end
            dLdu_s=[];
            for i=1:m
              e=zeros(1,m);e(i)=ei;
              dLdu_s=[dLdu_s (L(x,xr,u+e,ur,p,t,vdat)-L(x,xr,u-e,ur,p,t,vdat))/(2*ei)]; 
            end



end

 
% Return sensitivity right-hand side

dzS =[dfdx_s zeros(n,1);dLdx_s 0]*zS+[zeros(n) dfdu_s;zeros(1,n) dLdu_s];


flag = 0;new_data = [];


%------------------------------------------------------------------
% Gradient of the terminal cost
%------------------------------------------------------------------
function dx = dEdz(z,data)


[nt,np,n,m,ng,nb,M]=deal(data.sizes{1:7});
nz=nt+np+M*n+(M-1)*m;                       % Length of the optimization variable 

e=data.options.perturbation.J;                                % pertubation size;


if nt; tf=z(1); else tf=data.tf; end
% Extract p
p=z(nt+1:nt+np);
x0=z(nt+np+1:nt+np+n);
u0=z(nt+np+n+1:nt+np+n+m);
xf=z(end-n+1:end);
uf=z(end-m-n+1:end-n);


[L,E]=data.functions{1:2};
vdat=data.data;
if data.FD.index.dE.flag
    x=reshape(data.map.Vx*z,n,M)';
    usp=data.map.Vu*z;
    u=reshape([usp;usp(end-m+1:end)],m,M)';
    
    [dL,dE]=gradCost(L,x(end,:),data.references.xr,u(end,:),data.references.ur,p',0,E,x0,xf,u0,uf,p,tf,data);
    dx=spalloc(1,nz,nt+np+2*(n+m));
    dx(data.costStruct.E==1)=[dE.dtf dE.dp dE.dx0 dE.du0 dE.duf dE.dxf];
else

            
% Compute contribution of E to the gradient

            dEdx_s=spalloc(1,nz,nt+np+2*(n+m));
             idx=data.FD.index.Ey;
             nfd=size(idx,2);                               

            etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
            ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
            exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

         for i=1:nfd
             dEdx_s(1,idx(i))=(E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),0,tf+etf(:,i),vdat)-...
             E(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),0,tf-etf(:,i),vdat))/(2*e);
         end
         dx=dEdx_s; 

end

%------------------------------------------------------------------
% Jacobian of the user-defined constraint function
%------------------------------------------------------------------

function [dgdp,dg] = dgdz(x,u,p,t,data)

[nt,np,n,m,ng]=deal(data.sizes{1:5});


ei=data.options.perturbation.J;
vdat=data.data;

            g=data.functions{4};
            
            dgdx=[];dgdu=[];dgdp=[];
              
            
            for i=1:n
                e=zeros(1,n);e(i)=ei;
                dgdx=[dgdx, ((g(x+e,u,p,t,vdat)-g(x-e,u,p,t,vdat))/(2*ei)).'];
            end

            for i=1:m
                e=zeros(1,m);e(i)=ei;
                dgdu=[dgdu ((g(x,u+e,p,t,vdat)-g(x,u-e,p,t,vdat))/(2*ei)).']; 
            end

           if np
              for i=1:np
               e=zeros(np,1);e(i)=ei;
               dgdp=[dgdp ((g(x,u,p+e,t,vdat)-g(x,u,p-e,t,vdat))/(2*ei)).'];
              end
           else
              dgdp=zeros(ng,np);  
           end  
           

dg=[dgdx dgdu];


%------------------------------------------------------------------
% Jacobian of the user-defined terminal constraint function
%------------------------------------------------------------------
function db = dbdz(x0,x,u0,u,p,t,data)
[nt,np,n,m,ng,nb,M]=deal(data.sizes{1:7});

ei=data.options.perturbation.J;
vdat=data.data;
nz=n*M+m*(M-1)+np+nt;


if data.FD.index.db.flag
   [f,g,avrc,b]=deal(data.functions{3:6});
   [df,dg,dbz]=jacConst(f,g,x',u',p',t,b,x0,x,u0,u,p,t,data.t0,data);
   db=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
   db(data.costStruct.B==1)=[dbz.dtf dbz.dp dbz.dx0 dbz.du0 dbz.duf dbz.dxf];
    
else

                      
           b=data.functions{6};
           db=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
           idx=data.FD.index.b;nfd=size(idx,2);
           etf=ei*data.FD.vector.b.etf;ep=ei*data.FD.vector.b.ep;
           ex0=ei*data.FD.vector.b.ex0;eu0=ei*data.FD.vector.b.eu0;
           exf=ei*data.FD.vector.b.exf;euf=ei*data.FD.vector.b.euf;

          for i=1:nfd
            bp=b(x0+ex0(:,i),x+exf(:,i),u0+eu0(:,i),u+euf(:,i),p+ep(:,i),0,t+etf,vdat);
            bm=b(x0-ex0(:,i),x-exf(:,i),u0-eu0(:,i),u-euf(:,i),p-ep(:,i),0,t-etf,vdat);
            db=db+sparse(1:nb,idx(:,i),(bp-bm)/(2*ei),nb,nz);
          end

end


