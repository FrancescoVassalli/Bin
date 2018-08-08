#ifndef Photon_h
#define Photon_h
#include "Scalar.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <queue>
class Photon // Parton and myParticle are currently turned off to reduce dependency
{
  public:
    Photon(){}
    Photon(float _pT){
      this->pT = Scalar(_pT);
    }
    Photon(float _pT,float _phi){
      this->pT = Scalar(_pT);
      this->phi= Scalar(_phi);
    }
    Photon(double _pT,double _phi){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
    }
    Photon(double _pT,double _phi,bool process){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      direct=process;
    }
    Photon(float _pT,float _phi, float _eta){
      this->pT = Scalar(_pT);
      this->phi= Scalar(_phi);
      this->eta= Scalar(_eta);
    }
    Photon(int position,float _pT,float _phi, float _eta){
      this->pT = Scalar(_pT);
      this->phi= Scalar(_phi);
      this->eta= Scalar(_eta);
      this->position=position;
    }
    Photon(double _pT,double _phi, double _eta){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      this->eta= Scalar((float)_eta);
    }
    Photon(double _pT,double _phi, double _eta, bool process){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      this->eta= Scalar((float)_eta);
      direct=process;
    }
    Photon(int position, float* eT, float* phi, float* eta, int* id,int SIZE){
      etCone=.3;
      this->position =position;
      this->pT= Scalar(eT[position]);
      this->phi = Scalar(phi[position]);
      this->eta = Scalar(eta[position]);
      findIsoEt(phi,eta,eT,id,SIZE);
    }
/*    Photon(double _pT,double _phi, double _eta, bool process, std::queue<myParticle> all){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      this->eta= Scalar((float)_eta);
      direct=process;
      findIsoEt(all);
    }
    Photon(int position,double _pT,double _phi, double _eta, bool process, queue<myParticle> all){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      this->eta= Scalar((float)_eta);
      direct=process;
      this->position=position;
      findIsoEt(all);
    }
    Photon(double _pT,double _phi, double _eta,queue<myParticle> all){
      this->pT = Scalar((float)_pT);
      this->phi= Scalar((float)_phi);
      this->eta= Scalar((float)_eta);
      findIsoEt(all);
    }*/
    //needs ID to find ISO 
    /*Photon(int SIZE, int position, float* eT, float* phi, float* eta, bool direct,float eTCone){
      etCone=eTCone;
      this->position=position;
      pT=Scalar(eT[position]);
      this->phi=Scalar(phi[position]);
      this->eta=Scalar(eta[position]);
      this->direct=direct;
      findIsoEt(phi,eta,eT,id,SIZE);
      }*/
    Photon(int SIZE, int* id, float* eT, float* phi, float* eta, bool direct,float eTCone,float eTCut){
      etCone=eTCone;
      findPosition(SIZE,id,eT,eTCut);
      pT=Scalar(eT[position]);
      this->phi=Scalar(phi[position]);
      this->eta=Scalar(eta[position]);
      this->direct=direct;
      findIsoEt(phi,eta,eT,id,SIZE);
    }
    Photon(TLorentzVector tlv){
      pT=(float)tlv.Pt();
      eta=(float)tlv.Eta();
      phi=(float)tlv.Phi();
    }
    Photon(TLorentzVector tlv,int position){
      pT=(float)tlv.Pt();
      eta=(float)tlv.Eta();
      phi=(float)tlv.Phi();
      this->position=position;
    }
    ~Photon(){}
    Scalar getpT()const {
      return pT;
    }
    void setpT(float _pT){
      this->pT = Scalar(_pT);
    }
    Scalar getphi()const {
      return phi;
    }
    void setphi(float _phi){
      this->phi = Scalar(_phi);
    }
    Scalar geteta()const{
      return eta;
    }
    void seteta(float _eta){
      this->eta = Scalar(_eta);
    }
    /*void setParton(Parton p){
      parton=p;
    }*/
    bool isDirect()const{
      return direct;
    }
    int getPosition()const{
      return position;
    }
    /*float findIsoEt(queue<myParticle> all){
      isoEt=0;
      while(!all.empty()){
        if (inCone(all.front().gety(),all.front().getphi()))
        {
          isoEt+=all.front().geteT();
        }
        all.pop();
      }
      return isoEt;
    }*/

    float getIsoEt()const{
      return isoEt;
    }

    std::pair<Scalar,Scalar> getAngle()const{
      std::pair<Scalar,Scalar> r;
      r.first=phi;
      r.second=eta;
      return r;
    }
    inline float deltaR(float geta, float gphi)const{
      return TMath::Power((TMath::Power(TMath::Abs(geta-eta.value),2)+TMath::Power(deltaPhi(gphi,phi.value),2)),.5);
    }
    inline float deltaR(Photon p)const{ 
      return TMath::Power((TMath::Power(TMath::Abs(p.geteta().value-eta.value),2)+TMath::Power(deltaPhi(p.getphi().value,phi.value),2)),.5);
    }

  private:
    Scalar pT;
    Scalar phi;
    Scalar eta;
    bool direct;
    float isoEt;
    //Parton parton;
    float etCone = 0.3;
    int position;
    inline bool inCone(float geta, float gphi) const
    {
      if( sqrt(TMath::Power(TMath::Abs(geta-eta.value),2)+TMath::Power(deltaPhi(gphi,phi.value),2)) < etCone )
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    inline bool inCone(float geta, float gphi, int id) const
    {
      return (TMath::Power((TMath::Power(TMath::Abs(geta-eta.value),2)+TMath::Power(deltaPhi(gphi,phi.value),2)),.5) < etCone );
    }
    float findIsoEt(float* phi, float* eta, float* eT, int* id, int SIZE){
      isoEt=0;
      for(int i=0;i<SIZE;i++){
        if (inCone(eta[i],phi[i],id[i]))
        {
          isoEt+=eT[i];
        }
      }
      isoEt-=pT.value; //take the photon out
      return isoEt;
    }
    inline bool isPhoton(int id)const {
      return id==22;
    }
    int findPosition(int SIZE, int* id, float* et,float eTCut){
      for (int i = 0; i < SIZE; ++i)
      {
        if (isPhoton(id[i])&&et[i]>eTCut)
        {
          position=i;
          return position;
        }
      }
      return -1;
    }
    inline float deltaPhi(float i1, float i2)const{
      float r = TMath::Abs(i1-i2);
      if (r>TMath::Pi())
      {
        r= 2*TMath::Pi()-r;
      }
      return r;
    }
};
#endif
