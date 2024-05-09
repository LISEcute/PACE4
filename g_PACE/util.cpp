#include <stdlib.h>
#include <QtMath>
#include <QString>
#include <time.h>
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))


double pow_int(double par, int power);
double pow2(double par);
char *GetNextSymbol(char *s);
char *GetNextDelimeter(char *s);
char *GetIntFromString(int &V, char *s);
char *GetDoubleFromString(double &V, char *s);
double mzsqrt(double X);
QString ElementName(int IZ);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow_int(double par, int power)
{
double S=1;

for(int i=1; i<=power; i++) S*=par;
return S;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double pow2(double par)
{
return pow_int(par, 2);
}

//---------------------------------------------------------------------------
char *GetNextSymbol(char *s)
{
char c;

while(true) {
        c=*s;
        if(c=='\0' || c=='\n' || c=='\r') {return NULL;}

        if( !(c==' ' || c=='\t' || c==',')) break;
        s++;
        }
return s;
}
//---------------------------------------------------------------------------
char *GetNextDelimeter(char *s)
{
char c;

while(true) {
        c=*s;
        if(c=='\0' || c=='\n' || c=='\r')   {return NULL;}
        if( c==' ' || c=='\t' || c==',') break;
        s++;
        }
return s;
}

//---------------------------------------------------------------------------
char *GetIntFromString(int &V, char *s)
{
s=GetNextSymbol(s); if(!s) return NULL;
V=atoi(s);
return GetNextDelimeter(s);
}
//---------------------------------------------------------------------------
char *GetDoubleFromString(double &V, char *s)
{
s=GetNextSymbol(s); if(!s) return NULL;
V=atof(s);
return GetNextDelimeter(s);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double mzsqrt(double X) {
      X=max(0.,X);
      return sqrt(X);
      }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
QString ElementName(int IZ)
{
const char* symb=
    "n H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZn"
    "GaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPm"
    "SmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaU "
    "NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtB0B1B2B3B4B5B6B7B8B9C0C1C2C3C"
    "4C5C6C7C8C9D0  ";

    if(IZ<0 || IZ>130) IZ=131;

    QString ret;
    ret += (char)symb[IZ*2];
    char c2 = (char)symb[IZ*2+1];
    if( c2 != ' ') ret += c2;

    return ret;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

#define gdIA   16807
#define gdIM   2147483647
//#define gdAM   (1.0/gdIM)
#define gdIQ   127773
#define gdIR   2836
#define gdMASK 123459876

static long idum =13L;
static double gdAM=1.0/gdIM;

void InitRandom();
double MyRandom(void);


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void InitRandom(void)
{
    time_t now = time(0);
    struct tm *timep = localtime(&now) ;
    //gettime(&timep);

    idum =  timep->tm_year + (60*(60*timep->tm_hour + timep->tm_min)+timep->tm_sec);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double MyRandom(void)
{
long k;
double ans;
idum ^= gdMASK;  //XORing with MASK allows use of zero and other
k= idum/gdIQ;    //simple bit patterns for idum.

idum=gdIA*(idum-k*gdIQ)-gdIR*k; // Compute idum=(IA*idum) % IM without overows by Schrage's method.
if ( idum < 0) idum += gdIM;
ans=gdAM*idum; 			  // Convert idum to a floating result.
idum ^= gdMASK; 			  //  Unmask before return.
return ans;
}

