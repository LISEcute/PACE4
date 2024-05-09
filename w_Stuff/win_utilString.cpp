//#pragma GCC diagnostic ignored "-Wstringop-truncation"
//#pragma GCC diagnostic ignored "-Wstringop-overflow="

//#include "L_Init/myextern.h"
#include <QString>
#include <QMessageBox>
#include <QtMath>
#include <QDir>

#if !defined(liseStrcpyOS_h)
#include "w_Stuff/liseStrcpyOS.h"
#endif

//#include "w_Image_viewer.h"

//extern char DecimalPoint;
//extern QString LiseppPathForData;   // just for win_utilString -- not used in the code

char DecimalPoint='.';

char* qstrcpyL (char *dest,   size_t destsize, const QString &qs);
char* qstrncpyL(char *dest,  size_t destsize, const QString &qs, int n);


void getShortFileName(const char *source, char *dest, int opt);
void getShortFileName(QString &source, char *dest, int opt=0);
QString getShortFileName(const QString &source, int opt=0);

//void change_FileLisePath(QString &filename);
void change_FilePath(const QString &filename, QString &pathForData);
QString GetPathFromFilename(QString &source);

int GetExtension(QString &filename, char *extension);
int QGetExtension(const QString &filename, QString &extension);
int EraseExtention(QString &fileName);

int ErasePathFromFileName(QString &filename, const char *path);
int ErasePathFromFileName(QString &filename, const QString &path);
//int charErasePathFromFileName(char *filename, const QString &path);

int AddPathToFileName(QString &filename, const char *path);
int AddPathToFileName(QString &filename, const QString &path);
//char *strcat_allocate(char* dst, char *src);
//char *strcpy_new(char* dst, const char *src);
char *str_lowcase(char *src);
char *lower_string(char *str);

QString numberKillZero(double v, char ch,  int k, int m=0); // only for 'f'
QString numberKillExpZero(double v, int k, bool strong=true);
QString numberKillZero(QString str);
int kill_space(char *buffer, bool kill_all=true);
void kill_zero(char *ss, char fill);
void kill_zeroPoint(char *ss, char fill);
void killExpZero(char *ss, bool strong);
int changePathToUniversal(char *path);
int changePathToUniversal(QString& path);

FILE *mfopen(const QString& filename, const char* operand);

/*int QStringNumberPrintf(int k);
int QStringNumberPrintf(double k, char format, int precision);
QString QStringNumber(int k);
QString QStringNumber(double k, char format, int precision);
*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

char* qstrcpyL(char *dest, size_t destsize, const QString &qs)
{
  //QByteArray ba = qs.toLocal8Bit();
  //  strcpyL(dest,  destsize, ba.data());
  strcpyL(dest,  destsize, qs.toStdString().c_str());
  return dest;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char* qstrncpyL(char *dest,  size_t destsize, const QString &qs, int n)
{
  QString qtemp = qs.left(n);
  return qstrcpyL(dest, destsize, qtemp);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void change_FilePath(const QString& filename, QString &pathForData)
{
  int i_start  = filename.lastIndexOf('/')+1;
  int i_length = filename.length();  // the same as size()

  pathForData = filename;
  if(i_start > 0) pathForData.remove(i_start,i_length-i_start);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString GetPathFromFilename(QString &source)
{
  QString path = QDir::fromNativeSeparators(source);
  int i_start  = path.lastIndexOf('/')+1;
  int i_length = path.length();
  if(i_start > 0) path.remove(i_start,i_length-i_start);
  return path;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int ErasePathFromFileName(QString &filename, const QString &path)
{
  int n1 = filename.size();
  if(n1<=0) return -3;

  filename.remove(path, Qt::CaseInsensitive);
  int n2 = filename.size();

  if(n1!=n2) while(filename.at(0)==QChar('/'))                   //  n1 & n2 it is necessary in the case of UNIX no cut
    filename.remove(0,1);

  return filename.size();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int ErasePathFromFileName(QString &filename, const char *path)  //  should be finally erased
{
  int n1 = filename.size();
  if(n1<=0) return -3;

  //should use QString QDir::fromNativeSeparators(const QString &pathName)

  int n=-1;
  //const char *path = (const char *) pathQ;
  int L=(int)strlen(path);

  int lens = L+2;
  char *pb=new char[lens];
  char *pf=new char[lens];
  strcpyL (pb,lens,path);
  strcatL(pb,lens,"/");
  strncpy(pf,(const char*)filename.toStdString().c_str(),L+1); pf[L+1]='\0';
  lower_string(pb);
  lower_string(pf);

  if(strcmp(pb,pf)==0)
    {
      filename.remove(0,L+1);
      n=L+1;
    }

  delete[] pf;
  delete[] pb;

  return n;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int AddPathToFileName(QString &filename, const char *path)
{
  QString p=path;
  return AddPathToFileName(filename, p);
}

//--------------------------------------------------
int AddPathToFileName(QString &filename, const QString &path)
{
  //should use QString QDir::fromNativeSeparators(const QString &pathName)

  int n1 = filename.size();
  if(n1<=0) return -3;

  if(filename.indexOf("\\",0) >= 0  )    // exclude classical Windows case
    filename.replace("\\", "/");


  int n=filename.indexOf(':');

  if(n<0)
    {
      if(filename.indexOf("//",0) < 0  )    // exclude "net" case
        {
          if(!path.endsWith("/")) filename.insert(0,"/");
          filename.insert(0,path);
        }
    }


  if(filename.indexOf("\\",0) >= 0  )    // exclude classical Windows case
    filename.replace("\\", "/");


  n=filename.indexOf('\n');
  if(n>0)filename.remove(n,1);

  return n;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int EraseExtention(QString &filename)
{
  int n1 = filename.size();
  if(n1<=0) return -3;

  int n=filename.lastIndexOf(':');
  //int k=filename.lastIndexOf('\\');
  int k=filename.lastIndexOf('/');
  int p=filename.lastIndexOf('.');
  int size=filename.length();

  if(p>0 && p>k && p>n) filename.remove(p,size);
  else size=0;

  return size;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//char *strcat_allocate(char* dst, char *src)
//{
//  size_t size_swap = strlen(dst)+strlen(src)+1;
//  char *dst_swap=new char[size_swap];

//  if(dst) {
//      strcpyL(dst_swap,size_swap,dst);
//      strcatL(dst_swap,size_swap,src);
//    }
//  else    strcpyL(dst_swap,size_swap,src);

//  delete dst;
//  return dst_swap;
//}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW=
//char *strcpy_new(char* dst, const char *src)
//{
//  if(!src) return nullptr;   ///s ometimes doesnt work :(


//  if(dst) delete dst;
//  size_t L = strlen(src)+1;
//  dst = new char[L];
//  strcpyL(dst,L,src);
//  return dst;
//}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//bool GetVersionArray(char *str, int *array)
//{
//  array[0]=array[1]=array[2]=array[3]=array[4]=0;
//  // 0 -  Generation   7
//  // 1 -  Global       5
//  // 2 -  Current      152
//  // 3 -  Index	   (Generation+Global)
//  // 4 -  Index	   (Generation+Global+Current)
//  char *t;

//  array[0] = atoi(str);  t = strchr(str,'.');    if(!t) return false;
//  array[1] = atoi(++t);  t = strchr(t,  '.');

//  if(t) {array[2] = atoi(++t);}

//  array[3] = array[0]*100  + array[1];
//  array[4] = array[3]*1000 + array[2];
//  return true;
//}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getShortFileName(QString &source, char *dest, int opt)
{
  getShortFileName(source.toStdString().c_str(), dest, opt);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getShortFileName(const char *source, char *dest, int opt)
{
  QString temp(source);                       // 0 - only name
  // 1 - with extention
  int  i_start=temp.lastIndexOf('/');
  //int  i_start=temp->lastIndexOf('\\');

  if(i_start > 0) temp.remove(0,i_start+1);

  if(opt!=1) {
      i_start=temp.lastIndexOf('.');
      if(i_start > 0) temp.remove(i_start,temp.length());
    }

  qstrcpy(dest,(const char*)temp.toStdString().c_str());
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString getShortFileName(const QString &source, int opt)
{
  QString temp(source);
  if(temp.indexOf("\\",0) >= 0  )    // exclude classical Windows case
    temp.replace("\\", "/");

  int  i_start=temp.lastIndexOf('/');
  if(i_start > 0) temp.remove(0,i_start+1);

  if(opt!=1) {
      i_start=temp.lastIndexOf('.');
      if(i_start > 0) temp.remove(i_start,temp.length());
    }

  return temp;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int GetExtension(QString &filename, char *extension)
{
  QString temp(filename);

  int i_start=filename.lastIndexOf('.');

  temp.remove(0,i_start+1);

  if(temp.length() < 30) qstrcpy(extension,(const char*)temp.toStdString().c_str());
  else i_start = 0;

  return i_start;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int QGetExtension(const QString &filename, QString &extension)
{
  extension = filename;
  int i_start= extension.lastIndexOf('.');
  if(i_start>0) extension.remove(0,i_start+1);
  else extension="";
  return i_start;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int changePathToUniversal(char *path)
{
  int k=0;

  char *p=path;

  while(*p!='\x0' && *p!='\n')   // some gluke for "p"
    {
      if(*p == '\\') { *p = '/'; k++;}
      p++;
    }

  return k;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int changePathToUniversal(QString& path)
{
  path.replace('\\','/');
  return 1;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *str_lowcase(char *src)
{
  char *p=src;

  while(*p!='\x0' && *p!='\n')   // some gluke for "p"
    {
      if(*p >= 'A' && *p <= 'Z') *p += char('a' - 'A');
      p++;
    }
  return src;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
char *lower_string(char *str)
{
  char *str_lower=str;
  int i=0;
  while (str[i])
    {
      str_lower[i] = tolower(str[i]);
      i++;
    }
  return str_lower;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int sprintf_at_end(int Limit, char *buffer, const char *format, ...)
{
  int k = (int)strlen(buffer);

  if(k>=Limit) return EOF;

  char *t = &buffer[k];

  if(!t) return EOF;

  va_list args;

  va_start( args, format );
  k = vsprintf(t, format, args);
  va_end (args);

  return k;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString numberKillZero(double v, char ch,  int k, int m) // k - number after point, m = meaning points
{
  bool e_flag = tolower(ch) =='e';

  if(tolower(ch)=='g')
    {
      double va = fabs(v);
      int kn = qMax(0,k-1);
      double vn = pow(10,-kn);
      if(va>1e4 || va < vn ) e_flag=true;
      else ch='f';
    }

  if(e_flag) return QString::number(v,ch,k);


  if(tolower(ch)!='f' || k<0)
    QMessageBox::information(0,"\'numberZero\" function", "Wrong parameters");
  if(m<=0)   m=10;
  if(m<k-1 ) m=k-1; // up to decimal point

  if(qFabs(v)<1e-10) return QString("0");


  QString str(QString::number(v,ch,k));

  int size=str.length();
  int order = log10(qFabs(v));

  while(k>0)
    {
      if(k+order > m) {
          k--;
          str.remove(--size,1);
        }
      else break;
    }


  return numberKillZero(str);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString numberKillZero(QString str)
{
  int size=str.length();

  while(str.lastIndexOf(DecimalPoint)>0)   // do not remove 0 without decimal point
    {
      int i=str.lastIndexOf('0');
      if(i == size-1) str.remove(i,1);
      else {
          i=str.lastIndexOf(DecimalPoint);
          if(i == size-1) str.remove(i,1);
          break;
        }

      size--;
    }

  return str;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void killExpZero(char *ss, bool strong)
{
  QString  str(ss);
  size_t len = strlen(ss);

  int n = strong ? 2 : 1;

  for(int i=0; i<n; i++)
    {
      if(!str.endsWith("e-0"))
        str.replace("e-0","e-", Qt::CaseInsensitive);
      if(!str.endsWith("e+0"))
        str.replace("e+0","e+", Qt::CaseInsensitive);
    }

  qstrcpyL(ss, len+1, str);  // v.16.5.6
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString numberKillExpZero(double v, int k, bool strong)
{
  QString str = QString::number(v,'e',k);

  int n = strong ? 2 : 1;

  for(int i=0; i<n; i++)
    {
      if(!str.endsWith("e-0")) str.replace("e-0","e-", Qt::CaseInsensitive);
      if(!str.endsWith("e+0")) str.replace("e+0","e+", Qt::CaseInsensitive);
    }
  return str;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void kill_zero(char *ss, char fill)
{
  char *s,*s1, *s2;
  char sm[10];
  int i,l,lE=0, flagE=0, flagP=0;

  l=(int)strlen(ss);		if(l==0)return;
  s=strchr(ss,0);	if(!s)return;

  s1=strchr(ss,'E');

  if(s1) flagE=1;
  else {
      s1=strchr(ss,'e');
      if(s1)flagE=1;
    }

  s2=strchr(ss,'%');      if(s2!=nullptr){flagP=1; *s2=fill; s--; l--;}

  if(flagE)
    {
      lE=(int)strlen(s1);   l-=lE;
      s-=lE;
      if( lE==4 && *(s1+2)=='0') { *(s1+2)=*(s1+3); *(s1+3)=fill;  lE=3;}
      strncpyL(sm,10,s1,lE);
      for(i=0; i<lE; i++) *(s1++)=fill;
    }

  //if( (strchr(ss,'.')) != nullptr)
  if( (strchr(ss,DecimalPoint)) != nullptr)
    {
      while( 1 )
        {
          s--;
          l--;
          if(l==0)return;

          if( *s=='0'  )  {                  *s=fill;  continue; }
          //            else            { if(*s=='.')  *(s--)=fill;  break;    }
          else            { if(*s==DecimalPoint)  *(s--)=fill;  break;    }
        };
      s++;
    }

  if(flagE)
    {
      int k=10 -(sm-s);
      strncpyL( s, k, sm,lE); s+=lE;
    }

  if(flagP) *s='%';
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void kill_zeroPoint(char *ss, char fill)
{
  char *s,*s1, *s2;
  char sm[10];
  int i,l,lE=0, flagE=0, flagP=0;

  l=(int)strlen(ss);		if(l==0)return;
  s=strchr(ss,0);	if(s==nullptr)return;

  s1=strchr(ss,'E');

  if(s1!=nullptr)flagE=1;
  else {
      s1=strchr(ss,'e');
      if(s1!=nullptr)flagE=1; }

  s2=strchr(ss,'%');      if(s2!=nullptr){flagP=1; *s2=fill; s--; l--;}

  if(flagE)
    {
      lE=(int)strlen(s1);   l-=lE;
      s-=lE;
      if( lE==4 && *(s1+2)=='0') { *(s1+2)=*(s1+3); *(s1+3)=fill;  lE=3;}
      strncpyL(sm,10,s1,lE);
      for(i=0; i<lE; i++) *(s1++)=fill;
    }

  if( (strchr(ss, DecimalPoint)) != nullptr)
    {
      while( 1 )
        {
          s--;
          l--;
          if(l==0)return;

          if( *s=='0'  )  { *s=fill;  continue; }
          else            { if(*s==DecimalPoint)/**(s--)=fill;*/  break;    }
        };
      s++;
    }

  if(flagE) { strncpy( s, sm,lE); s+=lE;}
  if(flagP) *s='%';
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int kill_space(char *buffer, bool kill_all)
{
  char *ptr, *curr;

  curr=ptr=buffer;

  // kill_all
  if(kill_all)
    while(1)
      {
        if(*ptr=='\x0' || *ptr=='\n'){ *curr='\x0'; return (int)strlen(buffer);}
        if(*ptr!=' ')*(curr++)=*(ptr++);
        else ptr++;
      }
  else  {
      int L=(int)strlen(buffer);
      while(L>0)
        {
          if(buffer[L-1]==' ')
            {
              buffer[L-1]='\x0';
              L--;
            }
          else  return L;
        }
    }
  return 0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
FILE *mfopen(const QString& filename, const char* operand)
{

#if !defined(__CYGWIN__) && !defined(_WIN32) && !defined(_WIN64)

  FILE *f =   fopen(filename.toStdString().c_str(),operand);

#else

  QString woperand(operand);
  FILE *f = _wfopen(filename.toStdWString().c_str(), woperand.toStdWString().c_str() );

#endif

  return f;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
size_t strcpyL(char *dst, size_t dstsize, const char *src)  // strcpy_s
{
  size_t len = strlen(src);

  if(dstsize>0 && len>0)
    {
      size_t bl = (len < dstsize-1 ? len : dstsize-1);
      ((char*)memcpy(dst, src, bl))[bl] = 0;
      len=bl;
    }
  else len=0;

  return len;  // chars copied
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
size_t strncpyL(char *dst, size_t dstsize, const char *src, size_t count) // strncpy_s
{
  size_t len = strlen(src);
  len = std::min(len,count);

  if(dstsize>0 && len >0)
    {
      size_t bl = (len < dstsize-1 ? len : dstsize-1);
      ((char*)memcpy(dst, src, bl))[bl] = 0;
      len=bl;
    }
  else len=0;

  return len;  // char copied
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
size_t strcatL(char *dst, size_t dstsize, const char *src)
{
  size_t lenS = strlen(src);
  size_t lenD = strlen(dst);

  if(dstsize>0 && lenS>0)
    {
      size_t left = dstsize-1-lenD;
      size_t bl = (lenS < left ? lenS : left);
      ((char*)memcpy(&dst[lenD], src, bl))[bl] = 0;
      lenS = bl;
    }
  else lenS=0;

  return lenS;  // char copied
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW  temporary
/*int QStringNumberPrintf(int k) // for PACE4 streams
{
QString qstr;
qstr.number(k);
int m=printf("%s",qstr.toStdString().c_str());
return m;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int QStringNumberPrintf(double k, char format, int precision)
{
QString qstr;
qstr.number(k, format, precision);
int m=printf("%s",qstr.toStdString().c_str());
return m;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString QStringNumber(int k)
{
static QString qstr;
qstr.number(k);
return qstr;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString QStringNumber(double k, char format, int precision)
{
QString qstr;
qstr.number(k, format, precision);
return qstr;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
