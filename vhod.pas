procedure WriteInFile(OutFileAddress : String;Xka:masnk);
{������ �室���� 䠩�� ��� ��⥣���}
const
  	str:array [1..36] of string = (
{str[1]:=}'1   ����� (1 - �ண���; 2 - ���襭�� �ࡨ��)',
{str[2]:=}' ��砫쭠� �� (TT)' ,
{str[3]:=}' 1     ��᫮ ��⭨���',

{str[4]:=}'         �������' ,
{str[5]:=}'��砫�� ������ �ண���� (TT)' ,
{str[6]:=}'������ ������ �ண���� (TT) ' ,
{str[7]:=}'��� �뤠� (ᥪ)' ,
{str[8]:=}'1E-3   �訡�� ����让 ����� (��)' ,

{str[9]:=}'���������' ,
{str[10]:=}'obsvyb.in  ���� � ������ﬨ',
{str[11]:=}'5           ��᫮ �������',
{str[12]:=}'2822.172315  2201.434565  5279.171158  ���न���� ���ࢠ�ਨ (��)',
{str[13]:=}'0.10         ��筮��� ���襭�� (�室�����) (��)',
{str[14]:=}'1E-5 1E-8    ��ਠ樨 ��砫��� ���न��� � ᪮��⥩ (�� � ��/c)',
{str[15]:=}'0.1        �����⥫�, �����騩 �室������',
{str[16]:=}'2            ��砫�� �᫮��� (1 - �� 䠩��; 2 - ��������� � �ணࠬ��)',
{str[17]:=}'8000        �ਡ����⥫�� ࠤ���-����� (��) �� ������ 1-�� �������',

{str[18]:=}'��������������',
{str[19]:=}'   10.  ����ﭭ� 蠣 ��⥣�஢���� (ᥪ; ��� ����. ��ࠬ���)',
{str[20]:=}'   19   ���冷� ��⥣��� (�� 7 �� 39 �१ 4)',
{str[21]:=}'   10   ��ࠬ��� ��⥣���',
{str[22]:=}' 1000   ���ࢠ� �஬������� �뤠� �� �࠭ (� 蠣�� ��⥣�஢����)',

{str[23]:=}'����������',
{str[24]:=}' 12  12   ��ମ���� �����⥭樠�� (NM)',
{str[25]:=}'    1   �㭠',
{str[26]:=}'    1   �����',
{str[27]:=}'    1   ���⮢�� �������� � �� ��䥪�',
{str[28]:=}'0 0 0   ����⨢���᪨� ��䥪�� (�_0, �_1, �_2)',
{str[29]:=}'    1   �ਫ���',
{str[30]:=}'    0   �⬮���',
{str[31]:=}'  100.  ���� ᣮ࠭�� (��)',

{str[32]:=}'�������',
{str[33]:=}'  500.   ���� (��)',
{str[34]:=}'  0.5   ���頤� �������� �祭�� (�^2)',
{str[35]:=}'   2.   ����-� �������� ᮯ�⨢�����',
{str[36]:=}'   2.   ����-� ��ࠦ����');
{_________________________________________________________________________}

var
  InFile, OutFile 		   : Text;
  i,j	                 	   : integer;
begin {pr2}
assign(outf,OutFileAddress);
rewrite(outf);

writeln(outf,str[1]);
writeln(outf,year1,'  ',month1,'  ',day1,'  ',hour1,'  ',min1,'  ',sec1:2:3,'  ',str[2]);
writeln(outf,str[3]);
writeln(outf,' ',xka[1],'  ',xka[2],'  ',xka[3]);
writeln(outf,' ',xka[4],'  ',xka[5],'  ',xka[6]);
writeln(outf);

writeln(outf,str[4]);
writeln(outf,year1,'  ',month1,'  ',day1,'  ',hour1,'  ',min1,'  ',sec1:2:3,'  ',str[5]);
writeln(outf,year2,'  ',month2,'  ',day2,'  ',hour2,'  ',min2,'  ',sec2:2:3,'  ',str[6]);
writeln(outf,step,' ',str[7]);
writeln(outf,str[8]);
writeln(outf);

for i:=9 to 17 do
  writeln(outf,str[i]);
writeln(outf);
for i:=18 to 22 do
  writeln(outf,str[i]);
writeln(outf);
for i:=23 to 31 do
  writeln(outf,str[i]);
writeln(outf);
for i:=32 to 36 do
  writeln(outf,str[i]);
close(outf);
end;{pr2}
{------------------------------------------------------------------------}




{---------����� ����஭��� �ந�������----------------------------------}
  For k:=1 to nk do
   begin
     if k < 4 then varc:=varc1 else varc:=varc2;

     for l:=1 to nk do
         xv[l]:=xu[l];
     xv[k]:=xv[k]+varc;

    {1}     WriteInfile(FileOut,xv);

    {2}      {����� ��⥣���--10����-}
              {swapvectors; }
              exec('iszm_puc.exe','');
              {swapvectors;
             {---------------------------}

    {3}     ReadInFile(FileIn2,cm);
     for l:=1 to nm do
         Am[l,k]:=(cm[l]-xc[l])/varc;
   end;
{-------------------------------------------------------------------------}

