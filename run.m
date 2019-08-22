addpath('./utils');
input_decimal=[6785, 3, 12, 129, 1232, 123, 23, 234, 21421, 12]
r=0
l=0
re=0 %false
m=6
K=size(input_decimal, 2)
p=6
EbN0=100
patches=2^r
params_in=[]

%calculate length of messages
    if (re==0)
         B_patch = m*(m+3)/2 + p - 1;
    else
         B_patch = m*(m+1)/2 + p - 1;
    end
B = patches*B_patch - sum(l(2:end));


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
%%

sumpropfound=0
sumtiming=0
encoder=Encoder(r,l,re,m,p,K,EbN0,input_bits, B, patches)
[Y, parity] = encoder.chirrup_encode
decoder=Decoder(Y,r,l,parity,re,m,p,K,patches,params_in)
[output_bits, timing_trial] = decoder.chirrup_decode(Y,r,l,parity,re,m,p,K)
propfound_trial = compare_bits(input_bits,output_bits);
sumpropfound = sumpropfound + propfound_trial
sumtiming = sumtiming + timing_trial