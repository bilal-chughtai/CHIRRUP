%This script is self contained, and prototypes the general root of unity
%chirp ideas outlined in the report


m=4;
k=2;
zeta=exp(2*pi*i/2^k)
P1=randP(m,k)
P2=randP(m,k)


y=gen_general_chirp_no_b(P1,zeta)%+gen_general_chirp_no_b(P2,zeta)



kbin=dec2bin(k);
iterations=k-1;
P_recovered = zeros(m);
iteration=1

while iteration <= iterations
    P_staging = zeros(m);
    
    for basis_index=1:m

        power=2^(iterations-iteration)
        e=zeros(m,1);
        e(basis_index,1)=1;
        probs=probe_ye(y,e,power,P_recovered, zeta);
        probs=abs(probs);
        [value, location]=sort(-probs);
        location=location(1)-1;
        rowpart = dec2bin(location,m)=='1'
        P_staging(basis_index,:)=P_staging(basis_index,:)+2^(iteration-1).*rowpart

    end
    P_recovered=P_recovered+P_staging
    iteration=iteration+1
    
end

P_recovered - P1
P_recovered - P2
    







function [prb] = probe_ye(y,e,power,P_recovered, zeta)

        % probe_ye  probe data vector with error vector
        %
        % [prb,yprod] = probe_ye(y,e)
        %
        % given data vector y and "error value" e,
        % compute Hadamard transform prb  of yprod = conj(y(a)) * y(a+e)


        % Sina Jafarpour 28/3/08

            Mpow = length(y);
            M = log2(Mpow);  % this should be an integer

            % index vector (starting at 0)
            a = (dec2bin((0:Mpow-1))=='1')';
            %a

            %peel_off
            for j = 1:Mpow
                peel_off(j,1) = zeta^(-a(:,j)'*P_recovered*a(:,j));
            end
            
          %  peel_off = zeta.^-(e'*P_recovered*a)

            % a plus e values
            apebin = bitxor(a,e);
            %apebin
            
            apedec=2.^(M-1:-1:0)*apebin;
            %apedec
            
            y=y.*peel_off
            
            
            % y(a+e)
            yape = y(apedec+1);
            %[y yape]

            % product
            yprod = conj(y).*yape;
            

            yprod = yprod.^power;
            
            %i^(e'*P*e)*i.^(2*e'*P_recovered*a)
            
           %expected_yprod = exp(i*pi/4)^(e'*P*e)*(-1).^(e'*P1*a)
            
%            yprod-expected_yprod.'
           

            % use fast Walsh-Hadamard transform routine
            prb = fhtnat(yprod);

        end


function rm = gen_general_chirp_no_b(P,zeta)
    M = size(P,1);
    % M=length(b);
    % constructs a general Read-Muller code as above, but with no b term, and replacing i with a general root of unity, e^2i*pi/k. k denotes the degree of root of unity
    rm = zeros(2^M,1);
    a = zeros(M,1);
    for q = 1:2^M
        sum1 = a'*P*a;
        
        %sum2 = b*a;
        rm(q) = zeta^sum1; %* (-1)^sum2;
        % next a
        for ix = M:-1:1 %interesting binary counter
            if a(ix)==1
                a(ix)=0;
            else
                a(ix)=1;
                break;
            end
        end
    end

end


function P = randP(m,k)
% convert a bitstring (as an integer or as a logical vector) to P and b
% * updated for Kerdock codes:  no longer zero diagonal.
%
% [P,b] = bstr2Pb(M,[]);                randomly assigned P and b
% [P,b] = bstr2Pb(M,[1 0 1 0 ...])      take bits from vector elements
% [P,b] = bstr2Pb(M,int);               take bits from int

% began life as a private function of genRM.m
% SJS 6/11/07
% SJS 28/3/08 rehashed for Kerdock codes:
%             diagonal no longer constrained to be all-zeros
%             total # free elements thus M(M+3)/2


    MM = ((m^2+m)/2);
    bstr = randi([0 2^(k-1)-1],MM,1);
    P = zeros(m,m);
    ixl = find(tril(ones(m),0));
    P(ixl)=bstr;
    P=P';
    P(ixl)=bstr;

end

function x=fhtnat(data)
        %------------------------------------------------------
        %1D Natural(Hadamard)ordered Fast Hadamard Transform
        %------------------------------------------------------
        % Author: Gylson Thomas
        % e-mail: gylson_thomas@yahoo.com
        % Asst. Professor, Electrical and Electronics Engineering Dept.
        % MES College of Engineering Kuttippuram,
        % Kerala, India, February 2005.
        % copyright 2007.
        % This function implements the 1D natural(Hadamard)ordered Fast Hadamard Transform,
            N = pow2(floor(log2(length(data))));
            x = data(1:N);
            k1=N; k2=1; k3=N/2;
            for i1=1:log2(N)
                L1=1;
                for i2=1:k2
                    for i3=1:k3
                        i=i3+L1-1; j=i+k3;
                        temp1= x(i); temp2 = x(j);
                        x(i) = temp1 + temp2;
                        x(j) = temp1 - temp2;
                    end
                        L1=L1+k1;
                end
                    k1 = k1/2;  k2 = k2*2;  k3 = k3/2;
            end
            x=inv(N)*x; %Delete this line for inverse transform
        end
