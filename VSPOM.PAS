{$I-}
unit vspom;
interface
const sizeofext=sizeof(extended);
      sizeofint=sizeof(integer);
function pelv(p:pointer;i:integer):extended;
function pelvi(p:pointer;i:integer):integer;
procedure elv(p:pointer;i:integer;a:extended);
procedure elvi(p:pointer;i:integer;a:integer);
function order(x:extended):extended;
function norder(x:extended):integer;
function askans:char;
function shift(p:pointer;n:integer):pointer;
procedure inpint(a:string;var b:integer);
procedure seekln(var f:text;n:integer);
procedure time_indicate;
function pelm(n:integer;p:pointer;i,j:integer):extended;
procedure elm(n:integer;p:pointer;i,j:integer;a:extended);
procedure umatr(n,kod1,kod2:integer;a,b,c:pointer);
procedure umatrn(n,kod1,kod2:integer;pa,pb,pc:pointer);
procedure mvec(n:integer;ap,bp,cp:pointer);
procedure copfile(name1,name2:string);

implementation
uses crt,dos;
function shift(p:pointer;n:integer):pointer;
{‘¤¢¨£ ¥â ãª § â¥«ì ­  ¯¥à¥¬¥­­ãî ­  n ¯®§¨æ¨© à §¬¥à  EXTENDED}
var pext:^extended;
begin pext:=p;
      inc(pext,n);
      shift:=pext
end;

function askans:char;
var ch:char;
begin
repeat
ch:=upcase(readkey)
until (ch='Y') or (ch='N');
writeln(' ',ch);
askans:=ch
end;

function order(x:extended):extended;
var absx,n:extended;
begin absx:=abs(x); n:=1;
  if absx=0 then
     begin
     writeln('®àï¤®ª ­ã«ï ­¥®¯à¥¤¥«¥­, ­ ¦¬¨â¥ <ENTER>');
     readln;
     halt end;
  if absx>1 then begin while absx>1 do
                       begin absx:=absx/10;
                             n:=n*10 end;
                 order:=n end;
  if absx<1 then while absx<1 do
                 begin absx:=absx*10;
                       n:=n/10 end;
  if x>0 then order:=n else order:=-n
end;
function norder(x:extended):integer;
var absx:extended;
    n:integer;
begin absx:=abs(x); n:=0;
  if absx=0 then
     begin
     writeln('®àï¤®ª ­ã«ï ­¥®¯à¥¤¥«¥­, ­ ¦¬¨â¥ <ENTER>');
     readln;
     halt end;

  if absx>1 then begin while absx>1 do
                       begin absx:=absx/10;
                             n:=n+1 end;
                 norder:=n end;
  if absx<1 then while absx<1 do
                 begin absx:=absx*10;
                       n:=n-1 end;
  norder:=n
end;

function pelv(p:pointer;i:integer):extended;
var pext:^extended;
begin pext:=p;
      inc(pext,i-1);
      pelv:=pext^;
end;
function pelvi(p:pointer;i:integer):integer;
var pext:^integer;
begin pext:=p;
      inc(pext,i-1);
      pelvi:=pext^;
end;

procedure elv(p:pointer;i:integer;a:extended);
var pext:^extended;
begin pext:=p;
   inc(pext,i-1);
   pext^:=a;
end;
procedure elvi(p:pointer;i:integer;a:integer);
var pext:^integer;
begin pext:=p;
   inc(pext,i-1);
   pext^:=a;
end;

procedure inpint(a:string;var b:integer);
begin writeln(a);
 repeat
 if inoutres<>0 then
   begin inoutres:=0;
   writeln('¥à¥¬¥­­ ï ¤®«¦­  ¡ëâì æ¥«®©. ®¢â®à¨â¥ ¢¢®¤.')
   end;
 readln(b)
 until inoutres=0
end;

procedure seekln(var f:text;n:integer);
{“áâ ­ ¢«¨¢ ¥â ãª § â¥«ì â¥ªáâ®¢®£® ä ©«  ­  n-ãî áâà®ªã}
{¥à¢ ï áâà®ª  ¨¬¥¥â ­®¬¥à 0.}
var i:integer;
begin reset(f);
for i:=0 to n-1 do readln(f)
end;

procedure time_indicate;
var hour,min,sec,sec100:word;
begin
gettime(hour,min,sec,sec100);
writeln('time:',hour:3,min:3,sec:3,'.',sec100:2)
end;

function pelm(n:integer;p:pointer;i,j:integer):extended;
var pext:^extended;
begin pext:=p;
      inc(pext,n*(i-1)+j-1);
      pelm:=pext^;
end;
procedure elm(n:integer;p:pointer;i,j:integer;a:extended);
var pext:^extended;
begin pext:=p;
   inc(pext,n*(i-1)+j-1);
   pext^:=a;
end;

procedure mvec(n:integer;ap,bp,cp:pointer);
{“Œ†…ˆ… Œ€’ˆ–› € ‚…Š’}
type arr=array[0..0] of extended;
var a:^arr absolute ap;
    b:^arr absolute bp;
    c:^arr absolute cp;
    s:extended;
    i,j,n_1:integer;
begin n_1:=n-1;
      for i:=0 to n_1 do
      begin s:=0;for j:=0 to n_1 do s:=s+a^[n*i+j]*b^[j];
      c^[i]:=s end
end;
procedure umatr(n,kod1,kod2:integer;a,b,c:pointer);
{“Œ†…ˆ… Œ€’ˆ–
 kod=0 - ‘’‚…’‘’‚“™€Ÿ Œ€’ˆ–€ …‡ ˆ‡Œ……ˆŸ,
 kod=1 - ‘’‚…’‘’‚“™€Ÿ Œ€’ˆ–€ ’€‘ˆ“…’‘Ÿ}
var i,j,k:integer;
    s:extended;
begin if (kod1<>1)and(kod2<>1) then
      for i:=1 to n do for j:=1 to n do
      begin s:=0;
      for k:=1 to n do s:=s+pelm(n,a,i,k)*pelm(n,b,k,j);elm(n,c,i,j,s) end;

      if (kod1=1)and(kod2=1) then
      for i:=1 to n do for j:=1 to n do
      begin s:=0;
      for k:=1 to n do s:=s+pelm(n,a,k,i)*pelm(n,b,j,k);elm(n,c,i,j,s) end;

      if (kod1=1)and(kod2<>1) then
      for i:=1 to n do for j:=1 to n do
      begin s:=0;
      for k:=1 to n do s:=s+pelm(n,a,k,i)*pelm(n,b,k,j);elm(n,c,i,j,s) end;

      if (kod1<>1)and(kod2=1) then
      for i:=1 to n do for j:=1 to n do
      begin s:=0;
      for k:=1 to n do s:=s+pelm(n,a,i,k)*pelm(n,b,j,k);elm(n,c,i,j,s) end
end;

procedure umatrn(n,kod1,kod2:integer;pa,pb,pc:pointer);
{“Œ†…ˆ… Œ€’ˆ–
 kod=0 - ‘’‚…’‘’‚“™€Ÿ Œ€’ˆ–€ …‡ ˆ‡Œ……ˆŸ,
 kod=1 - ‘’‚…’‘’‚“™€Ÿ Œ€’ˆ–€ ’€‘ˆ“…’‘Ÿ}
type arr=array[0..0] of extended;
var a:^arr absolute pa;
    b:^arr absolute pb;
    c:^arr absolute pc;
    i,j,k,n_1:integer;
    s:extended;
begin n_1:=n-1;
      if (kod1<>1)and(kod2<>1) then
      begin
        for i:=0 to n_1 do for j:=0 to n_1 do
        begin s:=0;
        for k:=0 to n_1 do s:=s+a^[n*i+k]*b^[n*k+j]; c^[n*i+j]:=s;
        end;
      exit;
      end;

      if (kod1=1)and(kod2=1) then
      begin
        for i:=0 to n_1 do for j:=0 to n_1 do
        begin s:=0;
        for k:=0 to n_1 do s:=s+a^[n*k+i]*b^[n*j+k]; c^[n*i+j]:=s;
        end;
      exit;
      end;

      if (kod1=1)and(kod2<>1) then
      begin
        for i:=0 to n_1 do for j:=0 to n_1 do
        begin s:=0;
        for k:=0 to n_1 do s:=s+a^[n*k+i]*b^[n*k+j]; c^[n*i+j]:=s;
        end;
      exit;
      end;

      if (kod1<>1)and(kod2=1) then
      for i:=0 to n_1 do for j:=0 to n_1 do
      begin s:=0;
      for k:=0 to n_1 do s:=s+a^[n*i+k]*b^[n*j+k]; c^[n*i+j]:=s end
end;

procedure copfile(name1,name2:string);
const dim=65520;
var inp,out:file;
    p:^byte;
    n,i,ost:longint;
begin getmem(p,dim);
assign(inp,name1); reset(inp,1);
assign(out,name2); rewrite(out,1);
n:=filesize(inp) div dim;
ost:=filesize(inp) mod dim;
for i:=1 to n do
  begin blockread(inp,p^,dim);
  blockwrite(out,p^,dim)
  end;
blockread(inp,p^,ost);
blockwrite(out,p^,ost);
close(inp); close(out);
freemem(p,dim);
end;

end.