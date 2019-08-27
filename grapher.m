r = 0
l = 0
re = 1
m = 6
p = 2
sigma = 0
mode = "rand"
trials = 5

output=[]
K=[]   
prop=1  
mvalues=[6,7,8]
size(mvalues,2)
pvalues=[5,6,7]

for a = 1:size(mvalues, 2)
    m=mvalues(a)
    
    for b=1:size(pvalues,2)
        p=pvalues(b)
        i=1

        while prop > 0.05
            i=i+1;
            K=[K, i];
            disp(["K=", i]);    
            [B, prop] = run(r, l, re, m, p, sigma, mode, i, trials);
            output = [output, prop];
        end
        
        filename=strcat("tests/B", num2str(B),'r', num2str(r),'l', num2str(l),'r', num2str(r),'m', num2str(m),'p', num2str(p), 'trials', num2str(trials))
        save(filename, "K", "prop");
    end
end



x=strcat("tests/B", num2str(B),'r', num2str(r),'l', num2str(l),'r', num2str(r),'m', num2str(m),'p', num2str(p), 'trials', num2str(trials))
x(1)
plot(K, output)
xlabel('K')
ylabel('prop found')