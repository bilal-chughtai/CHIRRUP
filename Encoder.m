classdef Encoder
    properties
        r
        l
        re
        m
        p
        K
        EbN0
        input_bits
        B
    end

    methods
        function self = Encoder(r,l,re,m,p,K,EbN0,input_bits, B)
            self.r = r;
            self.l = l;
            self.re = re;
            self.m = m;
            self.p = p;
            self.K = K;
            self.EbN0 = EbN0;
            self.input_bits = input_bits;
            self.B=B
        end



        function [Y input_bits parity] = chirrup_encode(self)

        %chirrup_encode  Generates K random messages and performs CHIRRUP encoding
        %
        % r          0, 1 or 2; 2^r patches
        % l          length-2^r vector of parity check digits
        %            recommend: r=0: l=0, r=1: l=[0 15], r=2: l = [0 10 10 15]
        % re         logical: false=complex chirps, true=real chirps
        % m          size of chirp code (2^m is length of codewords in each slot)
        %            recommend m = 6, 7 or 8
        % p          2^p slots
        %            require p<=(m-r)(m-r+3)/2-1 for complex
        %                    p<=(m-r)(m-r+1)/2-1 for real
        % K          number of messages
        % EbN0       energy-per-bit (Eb/N0)
        %
        % Y            Y{p} is a 2^m x 2^p matrix of measurements for patch p
        % input_bits   B x K matrix of the K B-bit messages
        % parity       parity check codes generated in the tree encoding
        %
        % No. of messages is B = 2^r*[(m-r-p)(m-r-p+3)/2+p-1]-sum(l)  for complex
        %                          B = 2^r*[(m-r-p)(m-r-p+1)/2+p-1)-sum(l) for real
        %
        % AJT (12/9/18)

            patches = 2^self.r;
            parity = [];



            %generate random messages
           % input_bits = rand(B,K)>0.5;

            %tree encoding
            if (patches>1)
                [patch_bits parity] = self.tree_encoder(self.input_bits,B_patch,patches,l);
                patch_bits = permute(patch_bits,[2 1 3]);
            else
                patch_bits = self.input_bits';
            end

            flag = false;

            %generate measurements for each patch
            for patch = 1:patches

                sigma = sqrt(patches*2^self.m/(self.B*self.EbN0));
                Y{patch} = self.sim_from_bits(sigma,patch_bits(:,:,patch));


            end
        end



        function [patch_bits parity] = tree_encoder(self,N,n)

        % tree_encoder  Encodes a message into patches using parity check digits
        %
        % input_bits      B x K matrix of the K B-bit messages
        % N               number of bits per patch
        % n               number of patches
        % l               length-n vector of number of parity check digits
        %                 recommend: n=1: l=0, n=2: l=[0 15], n=4: l = [0 10 10 15]
        %
        % patch_bits       N x K x n tensor of the K N-bit messages in each patch
        %
        % Regardless of input, zero parity check digits are used for the first
        % patch, and a number of parity check digits is used for the last patch to
        % agree with the total number of bits B. Ensure that B + sum(l) = N x n.
        %
        % Code is based on 'A Coupled Compressive Sensing Scheme for Unsourced
        % Multiple Access' by Amalladinne et al. 2018 (arXiv 1806.00138)

            self.l(1) = 0;
            %l(n) = N*n - size(input_bits,1) - sum(l) + l(n);

            patch_bits(:,:,1) = input_bits(1:N,:);
            count = N;

            for i = 2:n

                patch_bits(1:N-l(i),:,i) = input_bits(count+1:count+N-l(i),:);
                count = count + N - l(i);
                parity{i} = double(rand(l(i),count)>0.5);
                patch_bits(N-l(i)+1:N,:,i) = mod(parity{i}*input_bits(1:count,:),2);

            end
        end





        function Y = sim_from_bits(self,sigma,bits)

        % sim_from_bits  Generates binary chirp measurements from bits
        %
        % re         logical: false=complex chirps, true=real chirps
        % m          size of chirp code (2^m is length of codewords in each slot)
        %            recommend m = 6, 7 or 8
        % K          number of messages
        % p          2^p slots
        %            require p<=(m-r)(m-r+3)/2-1 for complex
        %                    p<=(m-r)(m-r+1)/2-1 for real
        % sigma      SD of noise: sigma = sqrt(patches*2^m/(B*EbN0))
        % bits       k x 2^m matrix of bits to encode
        %
        % Y          length-2^m vector of measurements
        %
        % AJT (12/9/18)


            Y = zeros(2^self.m,2^self.p);

            for k = 1:self.K

                %compute slots
                comps = self.compute_slots(bits(k,:));

                %make (P,b) for each slot
                bits1 = [0 bits(k,:)];
                bits2 = [1 bits(k,:)];
                [Pee1,bee1] = makePb(self.re,self.m,bits1);
                [Pee2,bee2] = makePb(self.re,self.m,bits2);

                %generate binary chirp vector for each slot
                rm1 = gen_chirp(Pee1,bee1);
                rm2 = gen_chirp(Pee2,bee2);

                %add onto measurement
                Y(:,comps(1)) = Y(:,comps(1))+rm1;
                Y(:,comps(2)) = Y(:,comps(2))+rm2;
            end

            %add noise (Gaussian for real, Complex Gaussian for complex)
            if (self.re==0)
                Y = Y + repmat(sigma*(randn(2^self.m,1)+1i*randn(2^self.m,1)),[1 2^self.p]);
            else
                Y = Y + repmat(sigma*randn(2^m,1),[1 2^self.p]);

            end
        end


        function comps = compute_slots(self, bits)
            if (self.re==0)
                nMuse = self.m*(self.m+1)/2;
            else
                nMuse = self.m*(self.m-1)/2;
            end

            disp(bits)
            comps(1) = outofbinary(bits(end-self.p+1:end))+1;

            trans = bits(nMuse+self.m-1:-1:nMuse+self.m-self.p);
            if outofbinary(trans)==0
                trans(1) = 1;
            end
            comps(2) = outofbinary(mod(bits(end-self.p+1:end)+trans,2))+1;
        end



    end
end