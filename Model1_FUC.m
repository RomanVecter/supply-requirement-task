function [cc,rf0,w2,w,Q,q,S,M,R,L]=Model1_FUC(c,a,k)
cc=c;
rf0=F_rf0_c(c,a,k);
w2=F_w2_c(c,a,k);
rf=0.01:0.01:1;
w=zeros(1,length(rf));
Q=zeros(1,length(rf));
q=zeros(1,length(rf));
S=zeros(1,length(rf));
M=zeros(1,length(rf));
R=zeros(1,length(rf));
L=zeros(1,length(rf));
for i=1:1:(length(rf))
    rdp=rf(i);
    t=rdp/(1+rdp);
    p=c*(1+rdp);
    if p>=1
        fprintf('���ڴ�������ʹ��c*(1+rf)>=1,��ʱc=%0.1f,rf0=%0.4f,w2=%0.4f,rf=%0.3f\n',c,rf0,w2,rdp);
    end
    if rdp<rf0
        %%�����
        w(i)=F_fan(H_inverse(t,a,k),a,k);
        Q(i)=F_fan_inverse(p,a,k);
        q(i)=H_inverse(t,a,k);
        %%����������
        S(i)=c*F_fan_inverse(p,a,k);
        %%M���ֲ�������
        fun1=@(x) F_fan(x,a,k)-p;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(p,a,k);%%��������
        fun2=@(x) F_fan(x,a,k)-F_fan(H_inverse(t,a,k),a,k);%%���ֺ���2
        a2=0; %%��������
        b2=H_inverse(t,a,k);%%��������
        M(i)=integral(fun1, a1, b1)-integral(fun2, a2, b2);
        %%R���ֲ�������
        R(i)=integral(fun2, a2, b2);
        %%L���ֲ�������
        fun3=@(x) F_fan(x,a,k);%%���ֺ���2
        a3=0; %%��������
        b3=F_fan_inverse(p,a,k);%%��������
        L(i)=rdp*(F_fan(H_inverse(t,a,k),a,k)*H_inverse(t,a,k)-c*F_fan_inverse(p,a,k))+integral(fun3, a3, b3);
    else
        w(i)=w2;
        Q(i)=w2/c*F_fan_inverse(w2,a,k);
        q(i)=F_fan_inverse(w2,a,k);
        %%����������
        S(i)=w2*F_fan_inverse(w2,a,k);
        %%M���ֲ�������
        fun1=@(x) F_distribution(x,a,k);  %%���ֺ���1
        a1=F_fan_inverse(w2,a,k); %%��������
        b1=w2/c*F_fan_inverse(w2,a,k);%%��������
        M(i)=(w2/c-1)*F_fan_inverse(w2,a,k)-integral(fun1, a1, b1);
        %%R���ֲ�������
        fun1=@(x) F_fan(x,a,k)-w2;  %%���ֺ���1
        a1=0; %%��������
        b1=F_fan_inverse(w2,a,k);%%��������
        R(i)=integral(fun1, a1, b1);
        %%L���ֲ�������
        fun3=@(x) F_fan(x,a,k);%%���ֺ���2
        a3=0; %%��������
        b3=w2/c*F_fan_inverse(w2,a,k);%%��������
        L(i)=integral(fun3, a3, b3);
    end
end

end