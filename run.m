input_decimal=[6785, 3, 12, 129, 1232, 123 ]
r=0
l=0
re=0 %false
m=6
K=2
p=6
EbN0=100

%calculate length of messages
    if (re==0)
         B_patch = m*(m+3)/2 + p - 1;
    else
         B_patch = m*(m+1)/2 + p - 1;
    end
B = 2^r*B_patch - sum(l(2:end));


%verify inputs are valid
for input=1:length(input_decimal)
    if input_decimal(input) >= 2^B
        error(['input_bits must be in ', B ,' bits.'])
    end
end

input_bits=[]
%convert to binary
for input=1:length(input_decimal)
    binaryinput = intobinary(input_decimal(input),B);
    input_bits = [input_bits binaryinput];
end


encoder=Encoder(r,l,re,m,p,K,EbN0,input_bits, B)
encoder.chirrup_encode