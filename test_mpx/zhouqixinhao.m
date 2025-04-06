N = 8;
t = 1:1:N;
x = 1:0.1:N;
for ind = 1:16
    subplot(4,4,ind)
    hold on
    stem(t, sin(2* pi * ind * t / N)); 
    plot(x, sin( 2* pi * ind * x / N)); 
    ylim([-1, 1]); title(ind);
end
