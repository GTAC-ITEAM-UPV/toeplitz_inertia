function [Gn,Hn] = generadores_toep_block(A,s)
%  A matriz toeplitz block simétrica
%  s numero de bloques
% los bloques son toeplitz simétricos
% guardamos en una array de celdas la primera columna de cada 
% bloque Toeplitz
[~,n]=size(A);
m=n/s;
Tf=cell(s,s);
for j=1:s
    for i=1:s
        Tf{i,j}=A(1+(i-1)*m,1+(j-1)*m:m*j);
    end
end
Tc=cell(s,s);
for j=1:s
    for i=1:s
        Tc{i,j}=A(1+(i-1)*m:m*i,1+(j-1)*m);
    end
end
Gn=zeros(n,2*s);
Hn=zeros(n,2*s);
for i=1:s
    Gn(1+(i-1)*m,i)=1;
    Hn(m+(i-1)*m,i+s)=1;
end
%poner primera columna de H, ultima de G (columna 2*s)
for i=1:s
    v=Tf{1,i}(:);
    Hn(1+m*(i-1):m*i,1)=[v(2:end);0];
    v=Tf{i,s}(:);
    v=[0;-v(end:-1:2)];
    Gn(1+m*(i-1):m*i,2*s)=v;
end

%Poner resto de columnas de Hn

for col=2:s
    for i=1:s
      vp=Tf{col,i}(:);
      vminus=Tc{col-1,i}(:);
      vp=[vp(2:end);0];
      vminus=vminus(end:-1:1);
      Hn(1+m*(i-1):m*i,col)=vp-vminus;
    end  
end

%Poner resto de columnas de Gn
for col=s-1:-1:1
    for i=1:s
        vp=Tc{i,col+1}(:);
        vminus=Tf{i,col}(:);
        vminus=[0;vminus(end:-1:2)];
        Gn(1+m*(i-1):m*i,col+s)=vp-vminus;
    end
end   


end

