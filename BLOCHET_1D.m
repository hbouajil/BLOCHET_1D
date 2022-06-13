                              %% %%-%%-%% PROBLEME DE BLOCHET 1D %%-%%-%% %%
clear all
clc
close all
%% PARAMETRES D'ENTREE :
n=100;                % Nb de descritisation suivant x.
H=2.189;              % Inclinaison des plans.
Delta_x=1/n;          % Pas suivant x : ∆x=1/N.

tic
%-------------Inesialisation de pression--------------
for i=1:n+1
    p(i)=0;
end
%-------------   Calcule de pression   ---------------
tol=1e-7;
nr=inf;
itr=0;
sum=0;
var=1;
alpha=1.8;

for i=2:n
        h(i)=H-(H-1)*(i-1)*Delta_x;
        hm(i)=H-(H-1)*(i-1-0.5)*Delta_x;
        hp(i)=H-(H-1)*(i-1+0.5)*Delta_x;
end

while var>tol
    for i=2:n

        % Calcule numérique de débit de gradient pression par differance finis :
        p(i)=[(hp(i).^3/(hm(i).^3+hp(i).^3))*p(i+1)+(hm(i).^3/(hm(i).^3+hp(i).^3))*p(i-1)-6*Delta_x*((hp(i)-hm(i))/(hm(i).^3+hp(i).^3))]*alpha+(1-alpha)*p(i);
        % Calcule numérique de débit de couette :
        Q_c(i)=h(i)/2; 
        % Calcule numérique de débit de poiseuille :
        Q_p(i)=((p(i+1)-p(i-1))*h(i).^3)/(24*Delta_x);
        % Sommation de pression :
        sum=sum+p(i);
      
    end
    var=abs((sum-nr)/sum);
    nr=sum;
    sum=0;
    itr=itr+1;   
end
toc 
sumFf=0.0;
    for i=2:n   

sumFf=sumFf+( 1/hp(i) + hp(i)/2*((p(i+1)-p(i)) )/(Delta_x) )*Delta_x;

    end


%% Expressions numériques :
disp(['les résultats numériques : ']);

    p_max=max(p)                                         % Pression non-dimensionelle 
    F=nr*Delta_x                                         % Portance non-dimensionelle
    Qc=0.5*h(i)                                        ; % Débit couette non-dimensionelle
    Qp=(1/24)*(h(i).^3)*abs(p(i+1)-p(i-1))*(1/Delta_x) ; % Débit poiseuille non-dimensionelle
    Q=Qc+Qp                                              % Débit non-dimensionelle
    Ff=sumFf                                             %Pf

%% Expressions analytiques :
disp(['les résultats analytiques : ']);

    x=0.68;
    p_anal=6*(((H-1)*(1-x)*x)/((H+1)*(H-(H-1)*x).^2))  % Pression non-dimensionelle 
    F_anal=6*((log(H)/(H-1)^2)-(2/(H*H-1)))            % Portance non-dimensionelle
    Q_anal=H/(H+1)                                     % Débit non-dimensionelle
    Pf_anal=((H-1)/2)*F_anal+log(H)/(H-1)              % Puissance de frottement non-dimensionelle 

    
%% SOUS-Programme pour n=100

n1=20;
H1=2.189;
Delta_x1=1/n1;

for i=1:n1+1
    p1(i)=0;
end

tol1=1e-7;
nr1=inf;
itr1=0;
sum1=0;
var1=1;
alpha1=1.8;

for i=2:n1
        h1(i)=H1-(H1-1)*(i-1)*Delta_x1;
        hm1(i)=H1-(H1-1)*(i-1-0.5)*Delta_x1;
        hp1(i)=H1-(H1-1)*(i-1+0.5)*Delta_x1;
end

while var1>tol1
    for i=2:n1

        p1(i)=[(hp1(i).^3/(hm1(i).^3+hp1(i).^3))*p1(i+1)+(hm1(i).^3/(hm1(i).^3+hp1(i).^3))*p1(i-1)-6*Delta_x1*((hp1(i)-hm1(i))/(hm1(i).^3+hp1(i).^3))]*alpha1+(1-alpha1)*p1(i);
        
        sum1=sum1+p1(i);
      
    end
    var1=abs((sum1-nr1)/sum1);
    nr1=sum1;
    sum1=0;
    itr1=itr1+1;   
end

%% SOUS-Programme pour H=6

n2=100;
H2=2.189;
Delta_x2=1/n2;


for i=1:n2+1
    p2(i)=0;
end

tol2=1e-7;
nr2=inf;
itr2=0;
sum2=0;
var2=1;
alpha2=1.8;

for i=2:n2
        h2(i)=H2-(H2-1)*(i-1)*Delta_x2;
        hm2(i)=H2-(H2-1)*(i-1-0.5)*Delta_x2;
        hp2(i)=H2-(H2-1)*(i-1+0.5)*Delta_x2;
end

while var2>tol2
    for i=2:n2


        p2(i)=[(hp2(i).^3/(hm2(i).^3+hp2(i).^3))*p2(i+1)+(hm2(i).^3/(hm2(i).^3+hp2(i).^3))*p2(i-1)-6*Delta_x2*((hp2(i)-hm2(i))/(hm2(i).^3+hp2(i).^3))]*alpha2+(1-alpha2)*p2(i);
        

        sum2=sum2+p2(i);
      
    end
    var2=abs((sum2-nr2)/sum2);
    nr2=sum2;
    sum2=0;
    itr2=itr2+1;   
end

%% Affichage des résultats :

% H=2.189
subplot(1,3,1)
x_bar=linspace(0,1,n+1);
plot(x_bar,p)
grid on
xlabel('x=x/L')
ylabel('p_a_d_i_m_e_n_s_t_i_o_n_e_l')

%axis([0 1 0 0.3])

% H=1.1
subplot(1,3,2)
x_bar1=linspace(0,1,n1+1);
plot(x_bar1,p1)
grid on
xlabel('x=x/L')
ylabel('p_a_d_i_m_e_n_s_t_i_o_n_e_l')
% H=6
subplot(1,3,3)
x_bar2=linspace(0,1,n2+1);
plot(x_bar2,p2)
grid on
xlabel('x=x/L')
ylabel('p_a_d_i_m_e_n_s_t_i_o_n_e_l')
figure
plot(x_bar1,p1,x_bar2,p2)
