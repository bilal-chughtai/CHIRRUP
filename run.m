function [] = run(r, l, re, m, p, sigma, mode, userinput, trials)
%master testing function
%mode="manual" or "rand"
%if manual then input is K - the number of messages
%if rand, then input is input_decimal, an array of decimal inputs
    
    patches=2^r;
    params_in=[];
    
    %calculate length of messages
    if (re==0)
         B_patch = m*(m+3)/2 + p - 1;
    else
         B_patch = m*(m+1)/2 + p - 1;
    end
    B = patches*B_patch - sum(l(2:end));
    sumpropfound=0;
    sumtiming=0;

    %if mode is random generate some random inputs
    if mode=="rand"
        K=userinput;
        for trial = 1:trials
            input_bits = rand(B,K) > 0.5;
            [propfound_trial, timing_trial] = chirrup_test(r,l,re,m,p,K,sigma,input_bits,B,patches,params_in);
            sumpropfound=sumpropfound+propfound_trial;
            sumtiming=sumtiming+timing_trial;
        end
    end

    if mode=="manual"
        input_bits=[]
        K=length(userinput)
        for input=1:length(userinput)
            if userinput(input) >= 2^B
                error(['input values must be in at most', B ,' bits.'])
            end
            %convert to binary
            binaryinput = intobinary(userinput(input),B);
            input_bits = [input_bits binaryinput];
        end
        for trial = 1:trials
                [propfound_trial, timing_trial] = chirrup_test(r,l,re,m,p,K,sigma,input_bits,B,patches, params_in);
                sumpropfound=sumpropfound+propfound_trial;
                sumtiming=sumtiming+timing_trial;
        end
    end

    propfound = sumpropfound/trials
    ave_time = sumtiming/trials

end



function [propfound_trial, timing_trial] = chirrup_test(r,l,re,m,p,K,sigma,input_bits,B,patches,params_in)
    encoder=Encoder(r,l,re,m,p,K,sigma,input_bits, B, patches);
    [Y, parity] = encoder.chirrup_encode;
    decoder=Decoder(Y,r,l,parity,re,m,p,K,patches,params_in);
    [output_bits, timing_trial] = decoder.chirrup_decode(Y,r,l,parity,re,m,p,K);
    propfound_trial = compare_bits(input_bits,output_bits);
end
