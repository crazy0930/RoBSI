function mi = MI(a,b,Dim)
%aΪ�ο�Ӱ��bΪ����Ӱ��dimΪֱ��ͼά����ͨ������Ϊ255,127,63��
[M,N] = size(a);

dim = Dim;
hab = zeros(dim+1,dim+1);
ha = zeros(1,dim+1);
hb = zeros(1,dim+1);
%��aӰ����й�һ��
if max(max(a))~=min(min(a))
    a = (a-min(min(a)))/(max(max(a))-min(min(a)));
else
    a = zeros(M,N);
end
%��bӰ����й�һ��
if max(max(b))-min(min(b))
    b = (b-min(min(b)))/(max(max(b))-min(min(b)));
else
    b = zeros(M,N);
end
%����ת��Ϊ�������ݣ���߼����ٶ�
a = double(int16(a*dim))+1;
b = double(int16(b*dim))+1;
% a = int16(a)+1;
% b = int16(b)+1;
%ͳ��ֱ��ͼ
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
index = find(hab~=0);%�õ�hab�в�Ϊ0�����
p = hab/hsum;%��һ��
Hab = sum(sum(-p(index).*log(p(index))));%����������

hsum = sum(sum(ha));
index = find(ha~=0);
p = ha/hsum;
Ha = sum(sum(-p(index).*log(p(index))));%����aӰ�����

hsum = sum(sum(hb));
index = find(hb~=0);
p = hb/hsum;
Hb = sum(sum(-p(index).*log(p(index))));%����bӰ�����

mi = Ha+Hb-Hab;%���㻥��Ϣ