% main program of the simulation of jamming transition in HN3 using
% Wang-Landau algorithm
% Need to use with other functions files
% For more information about this model:   http://arxiv.org/abs/1409.8313

% Parameters
clear all
k=7;                    % number of levels
N=2^k+1;                % total number of sites
hn3=zeros(N,1);         % sites 
neighbors=zeros(N,3);   % their corresponding neighbors
f=exp(1);               % modifaction factor
logf=log(f);
minf=1+1e-8;
minlogf=log(minf);
maxMC=1e8;                  % max MC steps;
flatness=0.80;
initialMC=k*N*N;
% initialization: 
%   find the neighbor of every site and store them
%   randomly place particles according to 1-neighbor constraint

fprintf('\n  ************ Begin HN3Jamming: N= %d, flatness= %4.2f ************ \n', N,flatness);
datetime=clock();
fprintf('\n                       %d/%d/%d,%d:%d \n',datetime(1),datetime(2),datetime(3),datetime(4),datetime(5));

fprintf('\n  ****** Initializing & Finding possilbe max occupations ****** \n');
for i=0:N-1
    neighbors(i+1,:)=neighbor(i,N);
end

Nmax=k; % max number of particles
logg=zeros(Nmax+1,1);
H=zeros(Nmax+1,1);
N1=0;
N2=N1;
Ns=[0:1:Nmax]';

for i=1:initialMC
    j=ceil(rand*N);
    if hn3(j)==0
        if AddorNot(j,hn3,neighbors)==0
            i=i-1;
        else
            N2=N1+1;
            if N2+1>length(logg)
                Ns=[Ns;N2];
                N1=N2;
                logg=[logg;1];
                H=[H;1];
                Nmax=length(logg)-1;
                hn3(j)=1;
            elseif rand<exp(logg(N1+1)-logg(N2+1))
                hn3(j)=1;
                N1=N2;
                logg(N1+1)=logg(N1+1)+logf;
                H(N1+1)=H(N1+1)+1;
            else
                logg(N1+1)=logg(N1+1)+logf;
                H(N1+1)=H(N1+1)+1;
            end
        end
        
    else
        N2=N1-1;
        if rand<exp(logg(N1+1)-logg(N2+1))
            hn3(j)=0;
            N1=N2;
            logg(N1+1)=logg(N1+1)+logf;
            H(N1+1)=H(N1+1)+1;
        else
            logg(N1+1)=logg(N1+1)+logf;
            H(N1+1)=H(N1+1)+1;
        end
    end
end
Nmax=max(Ns);
if Nmax<73
    Nmax=73;
end
fprintf('\n    Max Occupation: %+5d; Min Occupation: %+5d \n', Nmax,0);
datetime=clock();
fprintf('                     %d/%d/%d,%d:%d \n',datetime(1),datetime(2),datetime(3),datetime(4),datetime(5));
           
fprintf('\n  ****** Random Sampling ****** \n');

logg=zeros(Nmax+1,1);
H=zeros(Nmax+1,1);
Ns=[0:1:Nmax]';
i=0;
% Wang-Landau algorithm 
while logf>minlogf
    i=i+1;
    j=ceil(rand*N);
    if hn3(j)==0
        if AddorNot(j,hn3,neighbors)==0
            i=i-1;
        else
            N2=N1+1;
            if N2+1>length(logg)
                Ns=[Ns;N2];
                logg=[logg;1];
                H=[H;1];
                N1=N2;
                hn3(j)=1;
                Nmax=length(logg)-1;
                fprintf('    Max number of particles extended to %d. \n',Nmax);
            elseif rand<exp(logg(N1+1)-logg(N2+1))
                hn3(j)=1;
                N1=N2;
                logg(N1+1)=logg(N1+1)+logf;
                H(N1+1)=H(N1+1)+1;
            else
                logg(N1+1)=logg(N1+1)+logf;
                H(N1+1)=H(N1+1)+1;
            end
        end
        
    else
        N2=N1-1;
        if rand<exp(logg(N1+1)-logg(N2+1))
            hn3(j)=0;
            N1=N2;
            logg(N1+1)=logg(N1+1)+logf;
            H(N1+1)=H(N1+1)+1;
        else
            logg(N1+1)=logg(N1+1)+logf;
            H(N1+1)=H(N1+1)+1;
        end
    end
    
    % check whether to change f
    if mod(i,N*k*5)==0
        if mean(H)/max(H)>flatness && min(H)/mean(H)>flatness
            logf=logf/2;
            f=sqrt(f);
            datetime=clock();
            fprintf('\n  %d/%d/%d,%d:%d \n',datetime(1),datetime(2),datetime(3),datetime(4),datetime(5));
            fprintf('  Monte Carlo steps finished: %d;\n', (i)/N);
            fprintf('  Modification factor f: %12.9f.\n',exp(logf));
            if logf<minlogf
                fprintf('  Modification factor is too small to make a difference! \n \n ');
                fprintf('\n  ************ End: N= %d, flatness= %5.2f ************ \n \n', N,flatness);
                break;
            end
            H=zeros(length(H),1);
        end     
    end
end
logg=logg-logg(1);

step=0.01;
mu=[4:step:12]';
Nmu=length(mu);
Z=zeros(Nmu,1);
nu=Z;
S=zeros(Nmu-1,1);
for i=1:Nmu
    Z(i)=sum(exp(logg).*exp(Ns*mu(i)));
    nu(i)=1/N*sum(Ns.*exp(logg).*exp(Ns*mu(i))/Z(i));
end
logZ=log(Z);
S=1/N.*logZ-mu.*nu;

subplot(321);
plot(Ns,logg)
xlabel('Number of occupied sites');
ylabel('log(g)');  
axis([min(Ns) max(Ns) min(logg) max(logg)*1.05]);

subplot(322);
plot(1./mu,nu);
xlabel('1/chemical potential');
ylabel('Packing fraction');
axis([min(1./mu)*0.99 max(1./mu)*1.01 min(nu)*0.99 max(nu)*1.01]);

subplot(323);
plot(mu,nu);
xlabel('chemical potential');
ylabel('Packing fraction');
axis([min(mu)*0.99 max(mu)*1.01 min(nu)*0.99 max(nu)*1.01]);

subplot(324);
plot(mu,S);
xlabel('chemical potential');
ylabel('Entropy density');
axis([min(mu)*0.95 max(mu)*1.01 0 max(S)*1.05]);


subplot(325);
plot(nu,S);
xlabel('Packing Fraction');
ylabel('Entropy density');
axis([min(nu)*0.995 max(nu)*1.005 min(S)*0.95 max(S)*1.05]);

