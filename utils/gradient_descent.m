function [ H,grad ] = gradient_descent(H,w,F,rho,P,A,Z)

[n,k] = size(H);
%H=reshape(x,n,k);
% dv:n*1
alpha = 0.01;    
num_iterations = 10;
%theta = H; 
theta = reshape(H,n*k,1);
%H_history = zeros(n*k, 1);
%H_ini = H;
for iter = 1:num_iterations
    Ht = H;
    dv=1./sqrt(max(sum(bsxfun(@times,H,w),2),eps));
    de=1./max(sum(H),eps);
    %H_Hini=H-(P'*A*Z)';
    %obj1=sum(sum(H_Hini.^2));
    DvH=bsxfun(@times,H,dv);
    YDvH=F'*DvH;
    DeW=w.*de;
    WDe2=DeW.*de;
    YDvHDeW=bsxfun(@times,YDvH,DeW);
    %obj4=rho*sum(sum(YDvHDeW.*YDvH));

    %obj=obj1-obj4;
    %obj_his(iter) = obj;

    grad1=2.*(H-(P'*A*Z)');
    YDvHDeW2=bsxfun(@times,YDvH,WDe2);
    grad3=repmat(sum(YDvH.*YDvHDeW2),n,1);
    WDeHDvY=YDvHDeW';
    HWDeHDvY=H*WDeHDvY;
    grad4=(sum(HWDeHDvY.*F,2).*(dv.^3))*w;
    grad5=-2.*bsxfun(@times,F*YDvHDeW,dv);
    grad=grad1+rho.*(grad3+grad4+grad5);
    grad=reshape(grad,n*k,1);
    theta = theta - alpha * grad; 
    H = reshape(theta,n,k);
    H(H<0) = 0;
    H(H>1) = 1;
    diffH = abs(norm(H - Ht,'fro')/norm(Ht,'fro'));

    if iter>2
    stop = diffH;
       if  stop <1e-1
           break
       end
    end
    %imagesc(H)
end

%plot(obj_his)
end