clc;
clear;
close all;

S=@(w)((1+4.5.*w.^2)./((-w.^4+3).^2+6.5.*w.^2));
wb=5;
w=-wb:.1:wb;
figure(1)
plot(w,S(w))
title('Targer PSD')

T=pi/wb;
num_terms=8; % for ARMA, select p=q=num_terms.
for ii=0:num_terms
    fun=@(w)((1+4.5.*w.^2)./((-w.^4+3).^2+6.5.*w.^2)).*cos(ii.*w*T);
    R(ii+1)=integral(fun,-wb,wb);
end
    
axis_p=0:num_terms;
figure(2)
scatter(axis_p,R) % an even function
title('Auto-correlation function')

B=-R(2:end)';
for ii=1:num_terms
    A(ii,ii:num_terms)=R(1:num_terms-ii+1);
end

for ii=2:num_terms
    for jj=1:(ii-1)
        A(ii,jj)=A(jj,ii);
    end
end

a=inv(A)*B;
G=sqrt((R(1)+R(2:end)*a)/2/wb); % b_0^hat

R_y_hatw=zeros(num_terms+1,1);
R_y_hatw(1)=2*wb*G;
for ii=1:num_terms
    for jj=1:ii
        R_y_hatw(ii+1)=R_y_hatw(ii+1)-a(jj)*R_y_hatw(ii-jj+1);
    end
end
figure(3)
scatter(axis_p,R_y_hatw)
title('R_y_w')

for ii=1:num_terms
    E(ii,ii:num_terms)=R_y_hatw(1:num_terms-ii+1);
end
E=-E;
C=E';
C=-C;
D=-2*wb*eye(num_terms);

MTX=[A,E;C,D];

VEC1=-R(2:end)';
VEC2=-R_y_hatw(2:end);
VEC=[VEC1;VEC2];

ab=inv(MTX)*VEC;

axis_p1=1:1:num_terms;
figure
scatter(axis_p1,ab(1:num_terms))

figure
scatter(axis_p1,ab(num_terms+1:end))


% plot the estimated PSD
syms w
kk=1:1:num_terms;
col_vec = exp(-1i.*kk.*w.*T);
fun1=abs(G+ab(num_terms+1:end)'*col_vec').^2./abs(1+ab(1:num_terms)'*col_vec').^2;%(G^2./abs(1+a'*col_vec').^2);

fun11=matlabFunction(fun1);
w_range=-wb:.05:wb;

figure
plot(w_range,S(w_range))
hold on
scatter(w_range,fun11(w_range))
title('ARMA simulation')
xlabel('\omega')
ylabel('PSD')