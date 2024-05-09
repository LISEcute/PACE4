#ifndef ParticleStatesH
#define ParticleStatesH

enum  EnumParticlesFlags {
 ep_n       = 0x0001,  //  neutrons
 ep_p       = 0x0002,  //  protons
 ep_a       = 0x0004,  //  alpha
 ep_g       = 0x0008,  //  gamma
 ep_maska   = 0x000F,  //   maska of particles
 ep_permit  = 0x0010,
 ep_opened  = 0x0020
 };

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

class ClassParticleFlags
{
public:

ClassParticleFlags(){ State = 0; UseGate=false; Zgate=10; Ngate=10;}

int GetParticleFlags() {return State;}
void SetParticleFlags(int t) {State=t;}


void ClearParticleState(unsigned int mask) {State &= (unsigned int)(~mask);}
void SetParticleState(unsigned int mask)   {State |= (unsigned int)(mask);}
bool IsParticleStateSet(unsigned int mask) { return (State & mask) ? true : false; }
bool CanBeOpenedParticleFile() {return IsParticleStateSet(ep_maska) && IsParticleStateSet(ep_permit);}
bool CanBeWirtten() {return IsParticleStateSet(ep_maska) && IsParticleStateSet(ep_opened);}
bool CheckGate(int Z, int N)
                {
                if(!UseGate) return true;
                if(Z!=Zgate) return false;
                if(N!=Ngate) return false;
                return true;
                }

void SetParticleStateValue(unsigned int mask, bool flag)
         {
         if(flag) SetParticleState(mask);
         else     ClearParticleState(mask);}


private:

int State;


public:
bool UseGate;
int Zgate,Ngate;
};

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
class  fusion_event {

public:

fusion_event(){Ex = x = y = z = particle_angle = particle_energy = 0; Z=N=A=J=MJ=0;}

  double Ex;
  int Z;
  int N;
  int A;
  int J;
  int MJ;
  double x;
  double y;
  double z;
  double particle_angle;
  double particle_energy;
//------------------------
  double Ad(){A=Z+N; return A;}
  int    Ai(){A=Z+N; return A;}
  bool testZN(fusion_event &c) { return (c.Z==Z && c.N==N);}
  bool testZN(fusion_event *c) { return (c->Z==Z && c->N==N);}
  void init  (fusion_event &c) {Z=c.Z; N=c.N; A=Z+N;}
  void init  (fusion_event *c) {Z=c->Z; N=c->N;A=Z+N;}
  void init  (int Zi, int Ni)  {Z=Zi; N=Ni;A=Z+N;}
  void init  (int Zi, int Ni, double Ei, int Ji, int MJi)  {Z=Zi; N=Ni; Ex=Ei; J=Ji; MJ=MJi;A=Z+N;}
//------------------------

};


#endif
