% check wether the empty site is allowed to add a particle
% this function is to serve HN3_WL_ToyModel.m
function flag = AddorNot(n,hn3,neighbors)
        flag=0;
        nonzero=nnz(neighbors(n,:));
        if n==2
            nonzero=nonzero+1;      % because of n==2, it's site 1; its neighbors are 0,2,3
        end
        
        totalNei=0;
        for i=1:nonzero
            totalNei=totalNei+hn3(neighbors(n,i)+1);
        end
        
        
        if totalNei==0                  % if no neighbor particle; we can add particle
           flag=1;
           return;
        
        elseif totalNei==1              % if there is 1 neighbor particle, check whether this neighbor particle has a neighbor particle
           for i=1:nonzero
               if hn3(neighbors(n,i)+1)==1
                en=neighbors(n,i);
               end
           end
           nonzero=nnz(neighbors(en+1,:));
           if en==1
               nonzero=nonzero+1;
           end
           
           totalNeien=0;
           for i=1:nonzero
               totalNeien=totalNeien+hn3(neighbors(en+1,i)+1);
           end
           
           if totalNeien==0
               flag=1;
               return;
           else
               flag=0;
               return;
           end   
        else
            flag=0;
            return;
        end

end
