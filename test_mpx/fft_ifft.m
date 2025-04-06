N = 20;
n = 1:N;
x = randn(N,1);

fx = fft(x);
y = ifft(fx);

% figure;hold on
% plot(x,'-k')
% plot(y,'ro')
% 
% cla;
% 
% stem(abs(fx))

for ind=n
    s(:,ind) = fx(ind)*exp(1j*2*pi*)
end
