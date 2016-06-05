% function main = main()
NMain = 1 ;
fid=fopen('GeneratedTestData.txt','wt');
for k =1 :NMain 
 N = 6; 
 M=rand(N);
 M=0.5*(M+M');
 L=100; %  magnitude
 for n=1:N
     for n1=1:N        
        M(n,n1)=L* M(n,n1);
        if ((n1>n+1) || (n1<n1-1))
            M(n,n1)=0;
            M(n1,n)=0;
        end
     end
 end

%  d = 200*rand(N,1); % The diagonal values
%  t = triu(bsxfun(@min,d,d.').*rand(N),1); % The upper trianglar random values
%  M = diag(d)+t+t.'; % Put them together in a symmetric matrix
 s = '';
 for i = 1:N
     s1 = '(';
     for j = 1:N
        s1 = strcat(s1, num2str(M(i,j)));
        if(j~=N)
            s1 = strcat(s1,',');
        end
     end
     s1 = strcat(s1,')');
     if(i~=N)
         s1 = strcat(s1,',');
     end
    
     s = strcat(s,s1);
 end
s = strcat('[',num2str(N),',',num2str(N),']','(',s,')','\n#');

fprintf(fid, s);
end
fclose(fid);
%  end
% function sign = randSign()
%     sign = randi(3,1);
%      if mod(sign,2) == 0
%         sign = 1;
%      else
%         sign  = -1;
%      end
%  end
%  
