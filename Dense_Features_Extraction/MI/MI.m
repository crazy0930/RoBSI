function mi = MI(a,b,Dim)
%a为参考影像，b为输入影像，dim为直方图维数，通常设置为255,127,63等
[M,N] = size(a);

dim = Dim;
hab = zeros(dim+1,dim+1);
ha = zeros(1,dim+1);
hb = zeros(1,dim+1);
%对a影像进行归一化
if max(max(a))~=min(min(a))
    a = (a-min(min(a)))/(max(max(a))-min(min(a)));
else
    a = zeros(M,N);
end
%对b影像进行归一化
if max(max(b))-min(min(b))
    b = (b-min(min(b)))/(max(max(b))-min(min(b)));
else
    b = zeros(M,N);
end
%重新转换为整型数据，提高计算速度
a = double(int16(a*dim))+1;
b = double(int16(b*dim))+1;
% a = int16(a)+1;
% b = int16(b)+1;
%统计直方图
for i=1:M
    for j=1:N
       indexx = a(i,j);
       indexy = b(i,j) ;
       hab(indexx,indexy) = hab(indexx,indexy)+1;
       ha(indexx) = ha(indexx)+1;
       hb(indexy) = hb(indexy)+1;
   end
end

hsum = sum(sum(hab));
index = find(hab~=0);%得到hab中不为0的序号
p = hab/hsum;%归一化
Hab = sum(sum(-p(index).*log(p(index))));%计算联合熵

hsum = sum(sum(ha));
index = find(ha~=0);
p = ha/hsum;
Ha = sum(sum(-p(index).*log(p(index))));%计算a影像的熵

hsum = sum(sum(hb));
index = find(hb~=0);
p = hb/hsum;
Hb = sum(sum(-p(index).*log(p(index))));%计算b影像的熵

mi = Ha+Hb-Hab;%计算互信息