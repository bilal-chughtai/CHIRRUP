r = 0
l = 0
re = 1
m = 6
p = 2
sigma = 0
mode = "rand"
trials = 30

output=[]
K=[]    

for i=1:20;
    K=[K, i];
    disp(["K=", i]);    
    [B, prop] = run(r, l, re, m, p, sigma, mode, i, trials)
    output = [output, prop];
end

x=strcat("B=", num2str(B),' patches=', num2str(2^r))
x(1)
plot(K, output)
xlabel('K')
ylabel('prop found')