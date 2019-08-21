classdef Encoder
    properties
        Y
        r
        parity
        re
        m
        p
        K
        params
        patches
    end

    methods
        function self = Decoder(Y,r,l,parity,re,m,p,K,params_in,patches)
            self.Y=Y
            self.r=r
            self.parity=parity
            self.re=re
            self.m=m
            self.p=p
            self.K=K
            self.patches=patches

            %default parameter values
            self.params.alpha = 0.1; %accept components if coefficient is within alpha of 1
            self.params.circuits = 5; %number of circuits in peeling decoder
            self.params.sparsity_factor = 3; %factor by which number of iterations exceeds
            %expected sparsity for a given slot
            self.params.tree_order = 3; %number of candidates per row in findPb tree search
            sumpropfound = 0 ;

            %read in parameter inputs
            if (nargin==9)
                if isfield(params_in,'alpha')
                    self.params.alpha = params_in.alpha;
                end
                if isfield(params_in,'circuits')
                    self.params.circuits = params_in.circuits;
                end
                if isfield(params_in,'sparsity_ratio')
                    self.params.sparsity_factor = params_in.sparsity_factor;
                end
                if isfield(params_in,'tree_order')
                    self.params.tree_order = params_in.tree_order;
                end
            end


        end




        function [output_decimal, timing] = chirrup_decode(self, Y,r,l,parity,re,m,p,K,params_in)

         %chirrup_decode  Performs CHIRRUP decoding
        %
        % Y          Y{p} is a 2^m x 2^p matrix of measurements for patch p
        % r          0, 1 or 2; 2^r patches
        % l          length-2^r vector of parity check digits
        %            recommend: r=0: l=0, r=1: l=[0 15], r=2: l = [0 10 10 15]
        % parity     parity check codes generated in the tree encoding
        %            leave empty if only 1 patch (r=0) is used
        % re         logical: false=complex chirps, true=real chirps
        % m          size of chirp code (2^m is length of codewords in each slot)
        %            recommend m = 6, 7 or 8
        % p          2^p slots
        %            require p<=(m-r)(m-r+3)/2-1 for complex
        %                    p<=(m-r)(m-r+1)/2-1 for real
        % K          number of messages
        % params     various; see below
        %
        % output_bits  B x K matrix of the K B-bit messages
        % timing       running time in seconds
        %
        % AJT (12/9/18)






            global outer_recov

            warning('off','MATLAB:rankDeficientMatrix');

            flag = false;

            for patch = 1:patches

                    tic
                    outer_recov = [];
                    count = 1;
                    Yp = Y{patch};
                    found = [];

                    %cycle through the slots repeatedly
                    for c = 1:self.params.circuits

                        for slot = 1:2^p

                            %run chirp reconstruction on given slot
                            sparsity_ratio = 3;
                            recov = self.chirp_rec(self,Yp(:,slot),ceil(params.sparsity_factor*self.K/2^(self.p-1)),slot,params);

                            [i j] = self.upper_indices();

                            if (~isempty(recov(1).P))

                                %find slot pairs for each recovered component
                                for r = 1:size(recov,2)
                                    Q=recov(self.r).P;
                                    for k = 1:length(i)
                                        Pvec(k) = Q(i(k),j(k));
                                    end
                                    Pbvec = [Pvec recov(self.r).b];
                                    already = 0;
                                    if (self.re==0)
                                        trans = Pbvec(self.m*(self.m+1)/2+self.m:-1:self.m*(self.m+1)/2+self.m-self.p+1);
                                    else
                                        trans = Pbvec(self.m*(self.m-1)/2+self.m:-1:self.m*(self.m-1)/2+self.m-self.p+1);
                                    end
                                    if outofbinary(trans)==0
                                        trans(1) = 1;
                                    end

                                    if (self.re==0)
                                        if recov(self.r).P(1,1)==1
                                            recov(self.r).P(1,1) = 0;
                                            recov(self.r).comps(2) = slot;
                                            recov(self.r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        else
                                            recov(self.r).comps(1) = slot;
                                            recov(self.r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        end
                                    else
                                        if recov(self.r).P(1,2)==1
                                            recov(self.r).P(1,2) = 0;
                                            recov(self.r).P(2,1) = 0;
                                            recov(self.r).comps(2) = slot;
                                            recov(self.r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        else
                                            recov(self.r).comps(1) = slot;
                                            recov(self.r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        end
                                    end

                                    %accept component if its coefficient is close to 1
                                    %if (abs(recov(r).c - 1)<params.alpha) %alternative condition
                                    if (real(recov(self.r).c)>1-params.alpha && real(recov(r).c)<1+params.alpha && abs(imag(recov(self.r).c))<params.alpha)
                                        for r2 = 1:count-1
                                            if (min(min(recov(r).P==outer_recov(r2).P)) && min(recovself.(self.r).b==outer_recov(r2).b))
                                                already = 1;
                                            end
                                        end
                                        if already==0
                                            outer_recov(count).P = recov(self.r).P;
                                            outer_recov(count).b = recov(self.r).b;
                                            outer_recov(count).c = recov(self.r).c;
                                            outer_recov(count).comps = recov(self.r).comps;
                                            count = count + 1;
                                        end
                                    end
                                end
                            end
                        end
                    end

                    %convert components back into bits
                    bit_matrix = find_bits(outer_recov,self.re,self.m,self.p);
                    if (size(bit_matrix,1)==0)
                        flag = 1;
                    end
                    processed_bits{patch} = bit_matrix';
                    timings(patch) = toc;

            end

                %patch together messages if necessary
                tic
                if (flag==1)
                    propfound_trial = 0;
                    sumpropfound = sumpropfound + propfound_trial;
                    timings(patches+1) = toc;
                    output_bits=[];
                else
                    if (patches>1)
                        output_bits = tree_decoder(processed_bits,l,self.parity);
                    else
                        output_bits = bit_matrix';
                    end
                    kmin = min(self.K,size(output_bits,2));
                    output_bits = output_bits(:,1:kmin);
                    timings(patches+1) = toc;
                end

                timing = sum(timings);
        end





        function recov = chirp_rec(self,y,nitlim,slot,params)

        % chirp_rec     Runs the chirp reconstruction algorithm, incorporating
        %               components found in other slots.
        %
        % y        measurement vector for given slot
        % re       logical: false=complex chirps, true=real chirps
        % nitlim   iteration limit
        % slot     slot number
        %
        % recov    multidimensional struct storing (P,b,c) for each component
        %          e.g. for component i: recov(1).P, recov(1).b, recov(1).c
        %
        % Sina Jafarpour (14/4/08) and AJT (12/9/18)

            global outer_recov

            foo.P = [];
            foo.b = [];
            foo.c = [];
            foo.comps = [];
            recov(1) = foo;
            allRM = [];
            ncomp = 0;

            %include components found in other slots
            count = 1;
            for r = 1:size(outer_recov,2)
                if outer_recov(r).comps(1)==slot
                    recov(count).P = outer_recov(r).P;
                    recov(count).b = outer_recov(r).b;
                    recov(count).c = outer_recov(r).c;
                    count = count + 1;
                elseif outer_recov(r).comps(2)==slot
                    recov(count).P = outer_recov(r).P;
                    if (re==0)
                        recov(count).P(1,1) = 1;
                    else
                        recov(count).P(1,2) = 1;
                        recov(count).P(2,1) = 1;
                    end
                    recov(count).b = outer_recov(r).b;
                    recov(count).c = outer_recov(r).c;
                    count = count + 1;
                end
            end

            %peel off components and recalculate residual
            if (count>1)
            for r = 1:size(recov,2)
                RM = gen_chirp(recov(r).P,recov(r).b);
                allRM = [allRM RM];
            end
            cr = allRM\y;
            y = y - allRM*cr;
            ncomp = size(recov,2);
            end

            M = log2(length(y));

            while (ncomp < nitlim && norm(y)>1e-3)

                [Phat bhat] = findPb(y,re,params);

                %determine component from P,b
                RM = gen_chirp(Phat,bhat);
                allRM = [allRM RM];

                %use all components & refine previous ests.
                cr = allRM\y;
                c = cr(end);
                allc = [recov(1:ncomp).c].';
                newc = [allc; 0]+cr;
                for q = 1:ncomp
                    recov(q).c = newc(q);
                end
                y = y - allRM*cr;

                %a component has been detected and handled.
                %record component parameters
                ncomp = ncomp+1;
                recov(ncomp).P = Phat;
                recov(ncomp).b = bhat;
                recov(ncomp).c = c;
                recov(ncomp).comp = RM;

            end
        end



        function [i j] = upper_indices(self)

            inds = 0:self.m^2-1;
            j = mod(inds,self.m) + 1;
            i = (inds - j + 1)/self.m + 1;
            if (self.re==0)
                upper = i<=j;
            else
                upper = i<j;
            end
            i = i(upper);
            j = j(upper);
        end

    end
end