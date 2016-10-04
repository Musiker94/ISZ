{$R-}
unit singul;
interface
procedure obrsing(kod,n:integer;const pa:array of extended;
                                var pobr,pu,plam:array of extended);
{
kod =0 - обычное обращение
    =1 - обощенное обращение
n - размер матрицы
pa -  входная матрица (поставлять a[1])
pu - матрица собственных векторов
plam - массив собственных чисел (по модулю)
pobr -обратная матрица (поставлять obr[1])
Пример:
var a,obr:array[1..6,1..6] of extended;
obrsing(0,6,a[1],obr[1]);
}
procedure sobstvectors(a,sigmap,sobst:pointer);
function pown(x:extended;n:integer):extended;
procedure polyn_mnk(kod,n,m:integer;const tp,pp:array of extended;
                                        var a:array of extended);
procedure svd(m,n:integer;const ap:array of extended;
                          var wp:array of extended;matu:boolean;
                          var up:array of extended;matv:boolean;
                          var vp:array of extended;
                          var ierr:integer;
                          var rv1p:array of extended);

procedure eigen(n,tmx:integer; pa,pt:pointer);
procedure eigen_sim(n:integer; pa,pd,pv:pointer);
{pa-входная матрица
результаты:
pd-массив собственных чисел
pv-матрица из собственных векторов
при вызове использовать @, т.е. eigen_sim(@a,@d,@v)}

function determ(n:integer;pa:pointer):extended;
procedure obrmatd(n:integer;pa,pb:pointer);

implementation
uses vspom;
procedure svd;
      label 190,210,270,290,360,390,410,430,460,475,490,510,520,540,
            565,575,580,650,1000,1001;
      type arr=array[1..1] of extended;
      var i,j,k,l,ii,kk,ll,i1,k1,l1,mn,its:integer;
          f,g,h,scale,anorm,s,c,x,y,z:extended;
          a:arr absolute ap;
          w:arr absolute wp;
          u:arr absolute up;
          v:arr absolute vp;
          rv1:arr absolute rv1p;
      function sign(a,b:extended):extended;
      begin if b>0 then sign:=abs(a);
            if b<0 then sign:=-abs(a);
            if b=0 then sign:=0;
      end;
    function dmax(a,b:extended):extended;
      begin
      if a>=b then dmax:=a
              else dmax:=b
      end;

      begin

      for i:=1 to m do for j:=1 to m do u[(i-1)*m+j]:=0;
      for i:=1 to n do for j:=1 to n do v[(i-1)*n+j]:=0;
      ierr:=0;
      for i:=1 to m do for j:=1 to n do u[(i-1)*m+j]:=a[(i-1)*n+j];
      g:=0;
      scale:=0;
      anorm:=0;
      for i:=1 to n do
      begin
      l:=i+1;
      rv1[i]:=scale*g;
      g:=0;
      s:=0;
      scale:=0;
      if i>m then goto 210;
      for k:=i to m do scale:=scale+abs(u[(k-1)*m+i]);
      if scale=0 then goto 210;
      for k:=i to m do
      begin
      u[(k-1)*m+i]:=u[(k-1)*m+i]/scale;
      s:=s+sqr(u[(k-1)*m+i])
      end;
      f:=u[(i-1)*m+i];
      g:=-sign(sqrt(s),f);
      h:=f*g-s;
      u[(i-1)*m+i]:=f-g;
      if i=n then goto 190;
      for j:=l to n do
      begin
      s:=0;
      for k:=i to m do s:=s+u[(k-1)*m+i]*u[(k-1)*m+j];
      f:=s/h;
      for k:=i to m do u[(k-1)*m+j]:=u[(k-1)*m+j]+f*u[(k-1)*m+i];
      end;
 190: for k:=i to m do u[(k-1)*m+i]:=scale*u[(k-1)*m+i];
 210: w[i]:=scale*g;
      g:=0;
      s:=0;
      scale:=0;
      if(i>m)or(i=n) then goto 290;
      for k:=l to n do scale:=scale+abs(u[(i-1)*m+k]);
      if scale=0 then goto 290;
      for k:=l to n do
      begin
      u[(i-1)*m+k]:=u[(i-1)*m+k]/scale;
      s:=s+sqr(u[(i-1)*m+k])
      end;
      f:=u[(i-1)*m+l];
      g:=-sign(sqrt(s),f);
      h:=f*g-s;
      u[(i-1)*m+l]:=f-g;
      for k:=l to n do rv1[k]:=u[(i-1)*m+k]/h;
      if i=m then goto 270;
      for j:=l to m do
      begin
      s:=0;
      for k:=l to n do s:=s+u[(j-1)*m+k]*u[(i-1)*m+k];
      for k:=l to n do u[(j-1)*m+k]:=u[(j-1)*m+k]+s*rv1[k]
      end;
 270: for k:=l to n do u[(i-1)*m+k]:=scale*u[(i-1)*m+k];
 290: anorm:=dmax(anorm,abs(w[i])+abs(rv1[i]));
      end;
      if  not matv then goto 410;
      for ii:=1 to n do
      begin
      i:=n+1-ii;
      if i=n then goto 390;
      if g=0 then goto 360;
      for j:=l to n do v[(j-1)*n+i]:=(u[(i-1)*m+j]/u[(i-1)*m+l])/g;
      for j:=l to n do
      begin
      s:=0;
      for k:=l to n do s:=s+u[(i-1)*m+k]*v[(k-1)*n+j];
      for k:=l to n do v[(k-1)*n+j]:=v[(k-1)*n+j]+s*v[(k-1)*n+i]
      end;
 360: for j:=l to n do
      begin
      v[(i-1)*n+j]:=0;
      v[(j-1)*n+i]:=0
      end;
 390: v[(i-1)*n+i]:=1;
      g:=rv1[i];
      l:=i
      end;
 410: if  not matu then goto 510;
      mn:=n;
      if m<n then mn:=m;
      for ii:=1 to mn do
      begin
      i:=mn+1-ii;
      l:=i+1;
      g:=w[i];
      if i=n then goto 430;
      for j:=l to n do u[(i-1)*m+j]:=0;
 430: if g=0 then goto 475;
      if i=mn then goto 460;
      for j:=l to n do
      begin
      s:=0;
      for k:=l to m do s:=s+u[(k-1)*m+i]*u[(k-1)*m+j];
      f:=(s/u[(i-1)*m+i])/g;
      for k:=i to m do u[(k-1)*m+j]:=u[(k-1)*m+j]+f*u[(k-1)*m+i]
      end;
 460: for j:=i to m do u[(j-1)*m+i]:=u[(j-1)*m+i]/g;
      goto 490;
 475: for j:=i to m do u[(j-1)*m+i]:=0;
 490: u[(i-1)*m+i]:=u[(i-1)*m+i]+1
      end;
 510: for kk:=1 to n do
      begin
      k1:=n-kk;
      k:=k1+1;
      its:=0;
 520: for ll:=1 to k do
      begin
      l1:=k-ll;
      l:=l1+1;
      if(abs(rv1[l])+anorm)=anorm then goto 565;
      if(abs(w[l1])+anorm)=anorm then goto 540
      end;
 540: c:=0;
      s:=1;
      for i:=l to k do
      begin
      f:=s*rv1[i];
      rv1[i]:=c*rv1[i];
      if(abs(f)+anorm)=anorm then goto 565;
      g:=w[i];
      h:=sqrt(f*f+g*g);
      w[i]:=h;
      c:=g/h;
      s:=-f/h;
      if  not matu  then continue;
      for j:=1 to m do
      begin
      y:=u[(j-1)*m+l1];
      z:=u[(j-1)*m+i];
      u[(j-1)*m+l1]:=y*c+z*s;
      u[(j-1)*m+i]:=-y*s+z*c
      end
      end;
 565: z:=w[k];
      if l=k then goto 650;
      if its=100 then goto 1000;
      its:=its+1;
      x:=w[l];
      y:=w[k1];
      g:=rv1[k1];
      h:=rv1[k];
      f:=((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
      g:=sqrt(f*f+1);
      f:=((x-z)*(x+z)+h*(y/(f+sign(g,f))-h))/x;
      c:=1;
      s:=1;
      for i1:=l to k1 do
      begin
      i:=i1+1;
      g:=rv1[i];
      y:=w[i];
      h:=s*g;
      g:=c*g;
      z:=sqrt(f*f+h*h);
      rv1[i1]:=z;
      c:=f/z;
      s:=h/z;
      f:=x*c+g*s;
      g:=-x*s+g*c;
      h:=y*s;
      y:=y*c;
      if  not matv then goto 575;
      for j:=1 to n do
      begin
      x:=v[(j-1)*n+i1];
      z:=v[(j-1)*n+i];
      v[(j-1)*n+i1]:=x*c+z*s;
      v[(j-1)*n+i]:=-x*s+z*c
      end;
 575: z:=sqrt(f*f+h*h);
      w[i1]:=z;
      if z=0 then goto 580;
      c:=f/z;
      s:=h/z;
 580: f:=c*g+s*y;
      x:=-s*g+c*y;
      if  not matu then continue;
      for j:=1 to m do
      begin
      y:=u[(j-1)*m+i1];
      z:=u[(j-1)*m+i];
      u[(j-1)*m+i1]:=y*c+z*s;
      u[(j-1)*m+i]:=-y*s+z*c
      end
      end;
      rv1[l]:=0;
      rv1[k]:=f;
      w[k]:=x;
      goto 520;
 650: if z>=0 then continue;
      w[k]:=-z;
      if  not matv then continue;
      for j:=1 to n do v[(j-1)*n+k]:=-v[(j-1)*n+k]
      end;
      goto 1001;
 1000:ierr:=k;
 1001:end;

procedure obrsing(kod,n:integer;const pa:array of extended;
                                var pobr,pu,plam:array of extended);
type arr=array[1..1] of extended;
var work:^arr;
    v,b,c:^arr;
    u:arr absolute pu;
    sigma:arr absolute plam;
    i,j,ierr,nom:integer;
    min:extended;
    a:arr absolute pa;
    obr:arr absolute pobr;
    out:text;
    st:pchar;
const vyksig=[1,2];
begin
       getmem(work,n*10);
       getmem(v,n*n*10);
       getmem(b,n*n*10);
       getmem(c,n*n*10);

{assign(out,'singul.out');
filesearch(st,'singul.out',' '); writeln('okkk');
if st[0]=#0 then rewrite(out) else append(out);}

{for i:=1 to n do
  begin for j:=1 to n do write(a[n*(i-1)+j]:11,' '); writeln end;}

svd(n,n,a,sigma,true,u,true,v^,ierr,work^);

{for i:=1 to n do write(out,sigma^[i]:20,' '); writeln(out);
close(out);}

for i:=1 to n do for j:=1 to n do b^[(i-1)*n+j]:=0;
if kod=0 then
  begin for i:=1 to n do b^[(i-1)*n+i]:=1/sigma[i];
  umatr(n,0,1,v,b,c);
  umatr(n,0,1,c,@u,@obr);
  end;
if kod=1 then
  begin min:=sigma[1];nom:=1;
  for i:=2 to n do if abs(sigma[i])<abs(min) then begin min:=sigma[i];
                                                           nom:=i end;
  for i:=1 to n do
      if (*i {in vyksig}<>nom*)
         sigma[i]>1.E-15
                            then b^[(i-1)*n+i]:=1/sigma[i]
                             else b^[(i-1)*n+i]:=0;
  umatr(n,0,1,v,b,c);
  umatr(n,0,1,c,@u,@obr);
end;
       freemem(work,n*10);
       freemem(v,n*n*10);
       freemem(b,n*n*10);
       freemem(c,n*n*10);
if (kod<>0)and(kod<>1) then
 begin writeln('kod должен быть 0 или 1,нажмите ENTER');
 readln;writeln('ха-ха');readln; halt end;
end;


procedure sobstvectors(a,sigmap,sobst:pointer);
type mat=array[1..6,1..6] of extended;
     mas6=array[1..6] of extended;
var u:mat;
    ierr:integer;
    work:mas6;
    out:text;
    i,j:integer;
    aa:^extended absolute a;
    ua:^extended absolute u;
    sobsta:^extended absolute sobst;
    sigma:^extended absolute sigmap;
begin
svd(6,6,aa^,sigma^,true,ua^,true,sobsta^,ierr,work);

{assign(out,'u_v.mat');rewrite(out);
for i:=1 to 6 do
 begin for j:=1 to 6 do write(out,u[i,j]:17,' ');writeln(out) end;
writeln(out);
for i:=1 to 6 do
 begin for j:=1 to 6 do write(out,sobst[i,j]:17,' ');writeln(out) end;
close(out)}
end;

function pown(x:extended;n:integer):extended;
var p:extended;
   an,i:integer;
begin  an:=abs(n);
p:=1;
for i:=1 to an do p:=p*x;
if n>=0 then pown:=p
        else pown:=1/p;
end;

procedure polyn_mnk(kod,n,m:integer;const tp,pp:array of extended;
                                        var a:array of extended);
type arr0=array[0..0] of extended;
     arr1=array[1..1] of extended;
var c,b,obr:^arr0;
    u,lam:^arr1;
    i,l,k,m1:integer;
    t:arr1 absolute tp;
    p:arr1 absolute pp;
    s:extended;
begin m1:=m+1;
      getmem(c,10*m1*m1);
      getmem(obr,10*m1*m1);
      getmem(b,10*m1);
      getmem(u,10*m1*m1);
      getmem(lam,10*m1);
for k:=0 to m do for l:=0 to k do
  begin s:=0;
  for i:=1 to n do s:=s+pown(t[i],k)*pown(t[i],l);
  c^[m1*k+l]:=s;
  end;
for k:=0 to m do
 begin
 for l:=k+1 to m do c^[m1*k+l]:=c^[m1*l+k];
 s:=0;
 for i:=1 to n do s:=s+pown(t[i],k)*p[i];
 b^[k]:=s;
 end;

obrsing(kod,m+1,c^,obr^,u^,lam^);

mvec(m+1,obr,b,@a);

      freemem(c,10*m1*m1);
      freemem(obr,10*m1*m1);
      freemem(b,10*m1);
      freemem(u,10*m1*m1);
      freemem(lam,10*m1);
end;

procedure eig_jacobi(n:integer;eivec:boolean; ap,dp,vp:pointer; var rot:integer);
type arr=array[1..1] of extended;
var sm,c,s,t,h,g,tau,theta,tresh:extended;
    p,q,i,j:integer;
    a:^arr absolute ap;
    d:^arr absolute dp;
    v:^arr absolute vp;
    b,z:^arr;

begin getmem(b,n*10);
      getmem(z,n*10);
if eivec then for p:=1 to n do for q:=1 to n do
              if p=q then v^[n*(p-1)+q]:=1 else v^[n*(p-1)+q]:=0;
for p:=1 to n do
  begin d^[p]:=a^[n*(p-1)+p]; b^[p]:=d^[p]; z^[p]:=0;
  end;
rot:=0;
for i:=1 to 50 do
  begin sm:=0;
  for p:=1 to n-1 do for q:=p+1 to n do sm:=sm+abs(a^[n*(p-1)+q]);
  if sm=0 then exit;
  if i<4 then tresh:=0.2*sm/sqr(n) else tresh:=0;
  for p:=1 to n-1 do for q:=p+1 to n do
    begin g:=100*abs(a^[n*(p-1)+q]);
    if(i>4)or(abs(d^[p])+g=abs(d^[p]))or(abs(d^[q])+g=abs(d^[q]))
    then a^[n*(p-1)+q]:=0
    else if abs(a^[n*(p-1)+q])>tresh then
      begin h:=d^[q]-d^[p];
      if abs(h)+g=abs(h) then t:=a^[n*(p-1)+q]/h else
        begin theta:=0.5*h/a^[n*(p-1)+q];
        t:=1/(abs(theta)+sqrt(1+sqr(theta)));
        if theta<0 then t:=-t;
        end;
      c:=1/sqrt(1+sqr(t));
      s:=t*c;
      tau:=s/(1+c);
      h:=t*a^[n*(p-1)+q];
      z^[p]:=z^[p]-h;
      z^[q]:=z^[q]+h;
      d^[p]:=d^[p]-h;
      d^[q]:=d^[q]+h;
      a^[n*(p-1)+q]:=0;
      for j:=1 to p-1 do
        begin g:=a^[n*(j-1)+p]; h:=a^[n*(j-1)+q];
        a^[n*(j-1)+p]:=g-s*(h+g*tau);
        a^[n*(j-1)+q]:=h+s*(g-h*tau);
        end;
      for j:=p+1 to q-1 do
        begin g:=a^[n*(p-1)+j]; h:=a^[n*(j-1)+q];
        a^[n*(p-1)+j]:=g-s*(h+g*tau);
        a^[n*(j-1)+p]:=h+s*(g-h*tau);
        end;
      for j:=q+1 to n do
        begin g:=a^[n*(p-1)+j]; h:=a^[n*(q-n)+j];
        a^[n*(p-1)+j]:=g-s*(h+g*tau);
        a^[n*(q-1)+j]:=h+s*(g-h*tau);
        end;
      if eivec then
      for j:=1 to n do
        begin g:=v^[n*(j-1)+p]; h:=v^[n*(j-1)+q];
        v^[n*(j-1)+p]:=g-s*(h+g*tau);
        v^[n*(j-1)+q]:=h+s*(g-h*tau);
        end;
      inc(rot);
      end;
    end;
  for p:=1 to n do
    begin b^[p]:=b^[p]+z^[p];
    d^[p]:=b^[p];
    z^[p]:=0;
    end;
  end;
freemem(b,n*10);
freemem(z,n*10);
end;

procedure eigen(n,tmx:integer; pa,pt:pointer);
type arr=array[1..1] of extended;
var a:^arr absolute pa;
    t:^arr absolute pt;
    eps,ep,aii,aij,aji,h,g,hj,aik,aki,aim,ami,tep,tem,d,c,e,akm,amk,
    cx,sx,cot2x,sig,cotx,cos2x,sin2x,te,tee,yh,den,tanhy,chy,shy,c1,c2,
    s1,s2,tki,tmi,tik,tim:extended;
    i,j,k,m,it,nless1:integer;
    mark,left,right:boolean;
label done,cont;
begin
mark:=false;
left:=false;
right:=false;
if tmx<>0 then
  begin {формирование единичной матрицы в массиве t}
  if tmx<0 then left:=true else right:=true;
  for i:=1 to n do
    begin t^[n*(i-1)+i]:=1;
    for j:=1+i to n do begin t^[n*(i-1)+j]:=0; t^[n*(j-1)+i]:=0 end;
    end;
  end;
ep:=5.4E-20; {машинное эпсилон}
eps:=sqrt(ep);
nless1:=n-1;
{Начало основного цикла. Возможно выполнение 50 итераций}
for it:=1 to 50 do
  begin
  if mark then begin tmx:=1-it; goto done end;
  {Проверка сходимости}
  for i:=1 to nless1 do
    begin aii:=a^[n*(i-1)+i];
    for j:=i+1 to n do
      begin aij:=a^[n*(i-1)+j]; aji:=a^[n*(j-1)+i];
      if (abs(aij+aji)>eps)or
         ((abs(aij+aji)>eps)and(abs(aii-a^[n*(j-1)+j])>eps)) then goto cont;
      end;
    end; {конец проверки сходимости}
  tmx:=it-1;
  goto done;
  {начало очередного преобразования}
  cont:
  mark:=true;
  for k:=1 to nless1 do for m:=k+1 to n do
    begin h:=0; g:=0; hj:=0; yh:=0;
    for i:=1 to n do
      begin aik:=a^[n*(i-1)+k]; aim:=a^[n*(i-1)+m];
      te:=sqr(aik); tee:=sqr(aim);
      yh:=yh+te-tee;
      if (i<>k)and(i<>m) then
        begin
        aki:=a^[n*(k-1)+i];
        ami:=a^[n*(m-1)+i];
        h:=h+aki*ami-aik*aim;
        tep:=te+sqr(ami);
        tem:=tee+sqr(aki);
        g:=g+tep+tem;
        hj:=hj-tep+tem;
        end;
      end;
    h:=h+h;
    d:=a^[n*(k-1)+k]-a^[n*(m-1)+m];
    akm:=a^[n*(k-1)+m];
    amk:=a^[n*(m-1)+k];
    c:=akm+amk;
    e:=akm-amk;
    if abs(c)<=ep then
      begin {матрица вращения R[i] - единичная}
      cx:=1;
      sx:=0;
      end else
      begin {вычисление элементов матрицы вращения R[i]}
      cot2x:=d/c;
      if cot2x<0 then sig:=-1 else sig:=1;
      cotx:=cot2x+(sig*sqrt(1+sqr(cot2x)));
      sx:=sig/sqrt(1+sqr(cotx));
      cx:=sx*cotx;
      end;
    if yh<0 then begin tem:=cx; cx:=sx; sx:=-tem end;
    cos2x:=sqr(cx)-sqr(sx);
    sin2x:=2*sx*cx;
    d:=d*cos2x+c*sin2x;
    h:=h*cos2x-hj*sin2x;
    den:=g+2*(sqr(e)+sqr(d));
    tanhy:=(e*d-h/2)/den;
    if abs(tanhy)<=ep then
      begin {матрица сдвига S[i] - единичная}
      chy:=1;
      shy:=0;
      end else
      begin {вычисление элементов матрицы сдвига S[i]}
      chy:=1/sqrt(1-sqr(tanhy));
      shy:=chy*tanhy;
      end;
    {вычисление элементов матрицы преобразования T[i]=S[i]*R[i]}
    c1:=chy*cx-shy*sx;
    c2:=chy*cx+shy*sx;
    s1:=chy*sx+shy*cx;
    s2:=-chy*sx+shy*cx;
    {проверка необходимости очередного преобразования}
    if (abs(s1)>ep)or(abs(s2)>ep) then
      begin {по крайней мере одно преобразование выполняется
            следующим образом}
      mark:=false;
      {левое преобразование}
      for i:=1 to n do
        begin aki:=a^[n*(k-1)+i]; ami:=a^[n*(m-1)+i];
        a^[n*(k-1)+i]:=c1*aki+s1*ami;
        a^[n*(m-1)+i]:=s2*aki+c2*ami;
        if left then
          begin {формирование левых собственных векторов}
          tki:=t^[n*(k-1)+i];
          tmi:=t^[n*(m-1)+i];
          t^[n*(k-1)+i]:=c1*tki+s1*tmi;
          t^[n*(m-1)+i]:=s2*tki+c2*tmi;
          end;
        end; {конец левого преобразования}
      {правое преобразование}
      for i:=1 to n do
        begin aik:=a^[n*(i-1)+k]; aim:=a^[n*(i-1)+m];
        a^[n*(i-1)+k]:=c2*aik-s2*aim;
        a^[n*(i-1)+m]:=-s1*aik+c1*aim;
        if right then
          begin {формирование правых собственных векторов}
          tik:=t^[n*(i-1)+k];
          tim:=t^[n*(i-1)+m];
          t^[n*(i-1)+k]:=c2*tik-s2*tim;
          t^[n*(i-1)+m]:=-s1*tik+c1*tim;
          end;
        end; {конец правого преобразования}
      end;
    end;
  end;
tmx:=50;
done:
end;

procedure eigen_sim(n:integer; pa,pd,pv:pointer);
type arr=array[1..1] of extended;
var a:^arr absolute pa;
    d:^arr absolute pd;
    b:^arr;
    i,j:integer;
begin
getmem(b,n*n*10);
for i:=1 to n do for j:=1 to n do b^[n*(i-1)+j]:=a^[n*(i-1)+j];

eigen(n,1,b,pv);

for i:=1 to n do d^[i]:=b^[n*(i-1)+i];
freemem(b,n*n*10);
end;

function determ(n:integer;pa:pointer):extended;
type mat=array[1..1] of extended;
label 2,3,5,8,14;
var a:^mat absolute pa;
    c:^mat;
    p,q,p1,q1:array[1..10] of integer;
    i,j,k,l,m,im,jm,e,s:integer;
    w:extended;

procedure emax(var im,jm:integer; pa:pointer);
var i,j,k:integer;
    am:extended;
    a:^mat absolute pa;
begin am:=0;
 for i:=1 to n do for j:=1 to n do
 if abs(a^[n*(i-1)+j])>abs(am) then
   begin am:=a^[n*(i-1)+j];im:=i;jm:=j end;
end;

begin getmem(c,n*n*10);
for i:=1 to n do for j:=1 to n do c^[n*(i-1)+j]:=0;
for l:=1 to n-1 do
  begin emax(im,jm,a);p[l]:=im;q[l]:=jm;
  i:=1;
5:k:=1;
14:if i=p[k] then goto 8;k:=k+1;
   if k<=l then goto 14;w:=-a^[n*(i-1)+jm]/a^[n*(im-1)+jm];
   for j:=1 to n do a^[n*(i-1)+j]:=a^[n*(i-1)+j]+w*a^[n*(im-1)+j];
 8:i:=i+1;if i<=n then goto 5;
  for j:=1 to n do begin c^[n*(im-1)+j]:=a^[n*(im-1)+j]; a^[n*(im-1)+j]:=0 end;
  end;
m:=0;
for i:=1 to n do begin for l:=1 to n-1 do if i=p[l] then goto 2
                       else if l=n-1 then p[n]:=i; 2:end;
for j:=1 to n do begin for l:=1 to n-1 do if j=q[l] then goto 3
                       else if l=n-1 then q[n]:=j; 3:end;
c^[n*(p[n]-1)+q[n]]:=a^[n*(p[n]-1)+q[n]];
for l:=1 to n do begin p1[l]:=p[l];q1[l]:=q[l] end;
m:=1;
for l:=1 to n do
  begin if p[l]<>l then
        begin m:=-m;e:=p[l];p[l]:=l;
         for j:=1 to n do for s:=l+1 to n do
         if (j=q[s]) and (p[s]=l) then p[s]:=e end;
        if q[l]<>l then
        begin m:=-m;e:=q[l];q[l]:=l;
         for i:=1 to n do for s:=l+1 to n do
         if (i=p[s]) and (q[s]=l) then q[s]:=e end;
  end;
w:=1;
for l:=1 to n do w:=w*c^[n*(p1[l]-1)+q1[l]];
determ:=m*w;
freemem(c,n*n*10);
end;

procedure obrmatd(n:integer;pa,pb:pointer);
type mat=array[1..1] of extended;
var a:^mat absolute pa;
    b:^mat absolute pb;
    c,vs:^mat;
    i,j,k,l,m,n1:integer;
    det,deter:extended;
begin n1:=n-1;
getmem(c,n1*n1*10);
getmem(vs,n*n*10);

for i:=1 to n do for j:=1 to n do vs^[n*(i-1)+j]:=a^[n*(i-1)+j];
det:=determ(n,vs);
 for i:=1 to n do for j:=1 to n do
 begin for k:=1 to i-1 do for l:=1 to j-1 do c^[n1*(k-1)+l]:=a^[n*(k-1)+l];
      for k:=i to n1  do for l:=j to n1  do c^[n1*(k-1)+l]:=a^[n*k+l+1];
      for k:=1 to i-1 do for l:=j to n1  do c^[n1*(k-1)+l]:=a^[n*(k-1)+l+1];
      for k:=i to n1  do for l:=1 to j-1 do c^[n1*(k-1)+l]:=a^[n*k+l];

  m:=1;
  for k:=1 to i+j do m:=-m;
  deter:=determ(n1,c);
  b^[n*(j-1)+i]:=m*deter/det;
 end;
freemem(c,n1*n1*10);
freemem(vs,n*n*10);
end;

end.
