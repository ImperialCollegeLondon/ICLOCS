function [ fzz ] = hessian_AN_F( df, Hf, fzz, M, n, nt, nz, f, X, U, P, t0, T, e, e2, DT, adjoint_f, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

      idx=data.FD.index.f; 
      nfd=size(idx,2);
      Tj=kron(T,ones(1,n));

    persistent solsave;  
    if isfield(vdat,'HesSave') && ~isfield(solsave,'ft') && ~isfield(solsave,'f')
        solsave=vdat.HesSave;
    end

    if isfield(solsave,'ft') && (size(solsave.ft,1)~=nfd || size(solsave.ft,2)~=nfd || size(solsave.ft{1},1)~=M)
        solsave = rmfield(solsave,'ft');
        if isfield(solsave,'ft_val_indicator')
            solsave = rmfield(solsave,'ft_val_indicator');
        end
    end
    if isfield(solsave,'f') && (size(solsave.f.fp1_save,1)~=nfd || size(solsave.f.fp1_save,2)~=nfd)
        solsave = rmfield(solsave,'f');
    end

    if data.FD.FcnTypes.Ftype==3 && (data.ProblemTypes.FixedTime || (~data.ProblemTypes.FixedTime && ~data.FD.FcnTypes.FTRelation)) && data.ProblemTypes.FixedTime
            if ~isfield(solsave,'ft')
                ft_save=cell(nfd,nfd);
                val_indicator=[];
                hessian_AN_F_saveFixTime;
                solsave.ft=ft_save;
                solsave.ft_val_indicator=val_indicator;
            end
            
            if ~isempty(df{1}) && isfield(data.options,'parfor') && data.options.parfor
                parfor k=1:size(solsave.ft_val_indicator,2)
                    i=solsave.ft_val_indicator(1,k);
                    j=solsave.ft_val_indicator(2,k);
                    ft=solsave.ft{i,j}.*adjoint_f;
                    fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
                end  
%                 C = cell(size(solsave.ft_val_indicator,2),1);
%                 parfor k=1:size(solsave.ft_val_indicator,2)
%                     i=solsave.ft_val_indicator(1,k);
%                     j=solsave.ft_val_indicator(2,k);
%                     ft=solsave.ft{i,j}.*adjoint_f;
%                     C{k}=[idx(:,i),idx(:,j),reshape(ft',M*n,1)];
%                 end  
%                 IJV = cell2mat( C );
%                 fzz=sparse(IJV(:,1),IJV(:,2),IJV(:,3),nz,nz);
            else
%                 for i=1:nfd
%                    for j=1:i
%                        if any(Hf{j,i},'all')
%                           ft=solsave.ft{i,j}.*adjoint_f;
%                           fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
%                        end
%                    end
%                 end  
                C = cell(nfd*nfd,1);
                k=1;
                for i=1:nfd
                   for j=1:i
                       if any(Hf{j,i},'all')
                          ft=solsave.ft{i,j}.*adjoint_f;
                          C{k}=[idx(:,i),idx(:,j),reshape(ft',M*n,1)];
                          k=k+1;
                       end
                   end
                end  
                IJV = cell2mat( C );
                fzz=sparse(IJV(:,1),IJV(:,2),IJV(:,3),nz,nz);
            end
    else
        hessian_AN_F_nominal;
    end

    
    global HesSave;
    if isempty(HesSave)
        HesSave=solsave;
    end


end

