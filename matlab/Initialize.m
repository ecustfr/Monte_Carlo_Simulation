function [r,e]=Initialize(n,box,Length,type)
    if type=='fcc' %为什么可以如此设置晶格
        %nc 个粒子一个晶格
        nc = round( (n/4)^(1/3) );
        cell = box/nc;
        box2 = box/2;
        r = zeros(n,3);
        e = zeros(n,3);
        r_fcc = [  0.25,0.25,0.25; 0.25,0.75,0.75; 0.75,0.75,0.25; 0.75,0.25,0.75];
        e_fcc = [1,1,1; 1,-1,-1; -1,1,-1; -1,-1,1]*sqrt(1/3); % 硬棒分子的取向，也就是方位角
        
        i=1;
        for ix=1:nc
            for iy=1:nc
                for iz=1:nc
                    for a =1:4
                        r(i,:) = r_fcc(a,:)+[ix-1,iy-1,iz-1];
                        r(i,:) = r(i,:)*cell;
                        r(i,:) = r(i,:)-box2;
                        e(i,:) = e_fcc(a,:);
                        
                        assert( ~overlap(  r(i,:),e(i,:),r((1:i-1),:),e((1:i-1),:),box ,Length  ), "Density too high")
                        i=i+1;
                    end
                end
            end
        end
        
    end
    
    if type=='ran' %随机分布
        %不太有用，除非这个相互作用力是soft的，或者密度很低
        %针对原子流体，只是length=0
        
        iter_max = 10000;
        
        r=zeros(n,3);
        
        e=zeros(n,3);
        
        
        disp('Random positions')
        
        for i=1:n
            iter=0;
            while true
               r(i,:)=(rand([1,3])-0.5)*box; 
               if ~overlap(r(i,:),e(i,:),r((1:i-1),:),e((1:i-1),:),box,Length)
                   break;
               end
               iter = iter+1;
               assert(iter<=iter_max,"Too many iteractions");
            end
        end
    end
end

function W=overlap(ri,ei,r,e,box,ell)
    %先不考虑硬球模型，ell==0
    
    W=false;
    tol=1e-6;
    %[nj,d]=size(r);
    if numel(r)==0
        return;
    end
    
    if ell<tol %对球状分子考虑 注：LJ并不是点状，只有球状或棒状
        rij = ri - r;
        rij = rij - round(rij/box)*box; %周期性边界条件
        rij_sq = sum(rij.^2 ,2);
        
        if sum(find(rij_sq<1))>0 %重叠判据有些奇怪
            W=true;
            return;
        end
    end
    
    %棒状分子被写在下面
    
end
        
        
        
        
        
    

    
    
    

    