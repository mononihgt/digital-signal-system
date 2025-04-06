clear;clc
N=10;
n=1:N;
x=randn(N,1);


fx=fft(x);
for ind=1:N
  s(:,ind)=fx(ind)*exp(1j*2*pi*(ind-1)*[n-1]/N)/N;
end

%%
figure;hold on
plot(x,'k')
plot(real(sum(s,2)),'ro')
plot(imag(sum(s,2)),'kx')

%%
figure;hold on
for ind=1:N
  plot(real(s(:,ind))+ind*0.5)
end

%%
figure;hold on
subplot 611
stem(s(:,1),'r','linewidth',2);axis off;ylim([-1 1])
subplot 612
stem(real(s(:,2)+s(:,10)),'r','linewidth',2);axis off;ylim([-1 1])
subplot 613
stem(real(s(:,3)+s(:,9)),'r','linewidth',2);axis off;ylim([-1 1])
subplot 614
stem(real(s(:,4)+s(:,8)),'r','linewidth',2);axis off;ylim([-1 1])
subplot 615
stem(real(s(:,5)+s(:,7)),'r','linewidth',2);axis off;ylim([-1 1])
subplot 616
stem(real(s(:,6)),'r','linewidth',2);axis off;ylim([-1 1])


figure;hold on
subplot 611
stem(imag(s(:,1)),'k','linewidth',2);axis off;ylim([-1 1])
subplot 612
stem(imag(s(:,2)+s(:,10)),'k','linewidth',2);axis off;ylim([-1 1])
subplot 613
stem(imag(s(:,3)+s(:,9)),'k','linewidth',2);axis off;ylim([-1 1])
subplot 614
stem(imag(s(:,4)+s(:,8)),'k','linewidth',2);axis off;ylim([-1 1])
subplot 615
stem(imag(s(:,5)+s(:,7)),'k','linewidth',2);axis off;ylim([-1 1])
subplot 616
stem(imag(s(:,6)),'k','linewidth',2);axis off;ylim([-1 1])

%%
s(:,4:8)=s(:,4:8)*0;
figure;hold on
plot(x,'k')
plot(real(sum(s,2)),'-ro')
plot(imag(sum(s,2)),'kx')