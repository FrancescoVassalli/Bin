#include "TLegend.h"
#include "TH1F.h"
#include <limits.h>
#include <queue>
#include <iostream>

namespace {
	short colors[7]={kRed,kBlue,kGreen+2,kMagenta+3,kOrange+4,kCyan+1,kMagenta-7};
	short styles[7]={kFullCircle,kOpenSquare,kFullTriangleUp,kFullDiamond,kFullCross,kFullStar,kOpenCircle};
	const char alphanum[] =
"0123456789"
"!@#$%^&*"
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz";
	}
	template<class T>
	T* zeroArray(int n,T* type){
		T* r = new T[n];
		for (int i = 0; i < n; ++i)
		{
			r[i] = 0;
		}
		return r;
	}
	float* zeroArray(int n){
		float* r = new float[n];
		for (int i = 0; i < n; ++i)
		{
			r[i] = 0;
		}
		return r;
	}
	template<class T>
	T* clone(T* a,int SIZE){
		T* r = new T[SIZE];
		for (int i = 0; i < SIZE; ++i)
		{
			r[i] = a[i];
		}
		return r;
	}
	/*template<class T>
	queue<T> clone(queue<T> *q){
		queue<T> r;
		while(!q->empty()){
			r.push(q->front());
		}

	}*/
	void makeBins(float* bins, int min, int nBins, float width){
		for(int i=0; i<=nBins;++i){
			bins[i] = min + width*i;
		}
	}
	float* makeBins(float min, float max, float width, int* nBins){ //nBins can't be NULL
		*nBins =(int)((max-min+.5)/width);
		float *bins= new float[*nBins];
		for(int i=0; i<= *nBins;++i){
			bins[i] = min + width*i;
		}
		return bins;
	}
	void makeMarkerNice(TH1F** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetMarkerStyle(styles[i-1]);
			(*h)->SetMarkerColor(colors[i-1]);
			h++;
		}
	}
	void makeMarkerNice(TGraph** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetMarkerStyle(styles[i-1]);
			(*h)->SetMarkerColor(colors[i-1]);
			h++;
		}
	}
	void makeMarkerNice(TGraphErrors** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetMarkerStyle(styles[i-1]);
			(*h)->SetMarkerColor(colors[i-1]);
			h++;
		}
	}
	void gNice(){
		gStyle->SetOptStat(0);
		gStyle->SetErrorX(0);
	}
	void fixOffset(TH1* plot){
		plot->GetYaxis()->SetTitleOffset(1);
		plot->GetXaxis()->SetTitleOffset(1);
	}
	void makeLineColors(TH1F** h, int n){
		for (int i = 1; i < n; ++i)
		{
			(*h)->SetLineColor(colors[i-1]);
			h++;
		}
	}
	void makeDifferent(TH1* h, int i){
		h->SetLineColor(colors[i]);
		h->SetMarkerStyle(styles[i]);
		h->SetMarkerColor(colors[i]);
	}
	void makeDifferent(TGraph* h, int i){
		h->SetLineColor(colors[i]);
		h->SetMarkerStyle(styles[i]);
		h->SetMarkerColor(colors[i]);
	}
	void makeDifferent(queue<TH1F*> hQ){
		int i=0;
		while(!hQ.empty()){
			hQ.front()->SetLineColor(colors[i]);
			hQ.front()->SetMarkerStyle(styles[i]);
			hQ.front()->SetMarkerColor(colors[i]);
			hQ.pop();
			i++;
		}
	}
	void makeDifferent(TH1F** h, int SIZE){
		for (int i = 0; i < SIZE; ++i)
		{
			h[i]->SetLineColor(colors[i]);
			h[i]->SetMarkerStyle(styles[i]);
			h[i]->SetMarkerColor(colors[i]);
		}
	}
	void makeLegendPoint(TLegend* tl, TH1F** h, int n, std::string *titles){
		for (int i = 0; i < n; ++i)
		{
			tl->AddEntry((*h++),titles++->c_str(),"p");
		}
	}
	void makeLegendLine(TLegend* tl, TH1F** h, int n, std::string *titles){
		for (int i = 0; i < n; ++i)
		{
			tl->AddEntry((*h++),titles++->c_str(),"l");
		}
	}
	void makeNiceHist(TH1* h){
		h->SetMarkerStyle(kFullCircle);
	}
	void axisTitles(TH1* h,std::string x, std::string y){
		h->GetYaxis()->SetTitle(y.c_str());
		h->GetXaxis()->SetTitle(x.c_str());
	}
	void axisTitles(TGraph* h,std::string x, std::string y){
		h->GetYaxis()->SetTitle(y.c_str());
		h->GetXaxis()->SetTitle(x.c_str());
	}
	void axisTitles(THStack* h,std::string x, std::string y){
		h->GetYaxis()->SetTitle(y.c_str());
		h->GetXaxis()->SetTitle(x.c_str());
	}
	void smallBorders(){
		gPad->SetBottomMargin(.1);
		gPad->SetTopMargin(.1);
	}
	void axisTitleSize(TH1F* h,float s){
		h->GetYaxis()->SetTitleSize(s);
		h->GetXaxis()->SetTitleSize(s);
	}
	void axisTitleSize(TH2F* h,float s){
		h->GetYaxis()->SetTitleSize(s);
		h->GetXaxis()->SetTitleSize(s);
	}
	void axisLabelSize(TH1F* h,float s){
		h->GetYaxis()->SetLabelSize(s);
		h->GetXaxis()->SetLabelSize(s);
	}
	template<class T>
	void axisTitleOffset(T* h, float s){
		h->GetYaxis()->SetTitleOffset(s);
		h->GetXaxis()->SetTitleOffset(s);
	}
	void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize,Double_t tsize)
	{
		//  Double_t tsize=0.032;
 		TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
 		marker->SetMarkerColor(color);  marker->SetNDC();
 		marker->SetMarkerStyle(mstyle);
 		marker->SetMarkerSize(msize);
 		marker->Draw();

 		TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
 		l.SetNDC();
 		l.DrawLatex(x,y,text);
	}
	void doubleZero(TGraph *h, float y, float x){
		h->GetYaxis()->SetRangeUser(0,y);
		h->GetXaxis()->SetRangeUser(0,x);
	}
	void doubleZero(TH1F *h, float y, float x){
		h->GetYaxis()->SetRangeUser(0,y);
		h->GetXaxis()->SetRangeUser(0,x);
	}
	void printHist(TH1F *h){
		cout<<h->GetName()<<":\n";
		for (int i = 0; i < h->GetNbinsX(); ++i)
		{
			cout<<"Bin"<<i<<":"<<h->GetBinContent(i)<<'\n';
		}
	}
	inline string getNextPlotName(int* plotcount){
		string r= "plot"+to_string(*plotcount);
		*plotcount = *plotcount+1;
		return r;
	}
	void normalizeTotal(TH1F** hlist,int SIZE){
		float sum=0;
		for (int i = 0; i < SIZE; ++i)
		{
			sum+=hlist[i]->Integral();
		}
		for (int i = 0; i < SIZE; ++i)
		{
			hlist[i]->Scale(1/sum);
		}
	}
	void normalizeBins(TH1F** hlist, int SIZE){ // note that this does not propagate the error
		for (int i = 1; i <= hlist[0]->GetNbinsX(); ++i)
		{
			double sum=0;
			for (int j = 0; j < SIZE; ++j)
			{
				sum+=hlist[j]->GetBinContent(i);
			}
			if(sum==0)continue;
			//cout<<"Sum"<<i<<":"<<sum<<'\n';
			for (int j = 0; j < SIZE; ++j)
			{
				hlist[j]->SetBinContent(i,hlist[j]->GetBinContent(i)/sum);
			}
		}
	}
	THStack* getStack(TH1F** hlist, int SIZE){
		string thisName = string(hlist[0]->GetName())+"Stack";
		THStack *r = new THStack(thisName.c_str(),hlist[0]->GetTitle());
		for (int i = 0; i < SIZE; ++i)
		{
			r->Add(hlist[i],"");
		}
		return r;
	}
	THStack* getStack(TH1F** hlist, int SIZE,string options){
		string thisName = string(hlist[0]->GetName())+"Stack";
		THStack *r = new THStack(thisName.c_str(),hlist[0]->GetTitle());
		for (int i = 0; i < SIZE; ++i)
		{
			r->Add(hlist[i],options.c_str());
		}
		return r;
	}
	THStack* getStack(TH1F** hlist, int SIZE,string options,string xTitle, string yTitle){ //must be drawn before it can be editted
		string thisName = string(hlist[0]->GetName())+"Stack";
		THStack *r = new THStack(thisName.c_str(),hlist[0]->GetTitle());
		for (int i = 0; i < SIZE; ++i)
		{
			r->Add(hlist[i],options.c_str());
		}
		return r;
	}
	

Scalar* yArray(queue<Point> q){
	Scalar *r = new Scalar[q.size()];
	int count =0;
	while (!q.empty())
	{
		r[count] = q.front().y;
		q.pop();
		count++;
	}
	return r;
}
Scalar* xArray(queue<Point> q){
	Scalar *r = new Scalar[q.size()];
	int count=0;
	while (!q.empty())
	{
		r[count] = q.front().x;
		q.pop();
		count++;
	}
	return r;
}
Scalar* xArray(Point* ps, int SIZE){
	Scalar *r = new Scalar[SIZE];
	for (unsigned int i = 0; i < SIZE; ++i)
	{
		r[i] = ps[i].x;
	}
	return r;
}
Scalar* yArray(Point* ps, int SIZE){
	Scalar *r = new Scalar[SIZE];
	for (unsigned int i = 0; i < SIZE; ++i)
	{
		r[i] = ps[i].y;
	}
	return r;
}
queue<Scalar> yQueue(queue<Point> in){
	queue<Scalar> out;
	while(!in.empty()){
		out.push(in.front().y);
		in.pop();
	}
	return out;
}
queue<Scalar> xQueue(queue<Point> in){
	queue<Scalar> out;
	while(!in.empty()){
		out.push(in.front().x);
		in.pop();
	}
	return out;
}
float* valueArray(Scalar* a, int SIZE){
	float *r = new float[SIZE];
	for (unsigned int i = 0; i < SIZE; ++i)
	{
		r[i] = a[i].value;
	}
	return r;
}
float* uncertaintyArray(Scalar* a, int SIZE){
	float *r = new float[SIZE];
	for (unsigned int i = 0; i < SIZE; ++i)
	{
		r[i] = a[i].uncertainty;
	}
	return r;
}
#ifndef Parton_h
#define Parton_h
#include "TLorentzVector.h"
class Parton
{
public:
	Parton(){}
	~Parton(){}
	Parton(int ID, float phi, float y){
		quark=isQuark(ID);
		this->phi=phi;
		this->y=y;
	}
	Parton(int ID, float phi, float eta, float eT, float e){
		quark=isQuark(ID);
		TLorentzVector tlv;
		tlv.SetPtEtaPhiE(eTTopT(eT,ID),eta,phi,e);
		this->phi=phi;
		y=tlv.Rapidity();
		this->eta = eta;
	}
	float getphi(){
		return phi;
	}
	float gety(){
		return y;
	}
	float geteta(){
		return eta;
	}
	bool getQuark(){
		return quark;
	}
protected:
	bool quark;
	float phi;
	float y;	
	float eta;
	bool isQuark(int ID){
		return TMath::Abs(ID)>0&&TMath::Abs(ID)<9;
	}
	float eTTopT(float eT, int ID){
		return TMath::Power(eT*eT-idToMass(ID)*idToMass(ID),.5);
	}
	float idToMass(int ID){
		ID=TMath::Abs(ID);
		switch(ID){
			case 2: return 0.0022;
				break;
			case 1: return 0.0047;
				break;
			case 3: return 0.096;
				break;
			case 4: return 1.28;
				break;
			case 5: return 4.18;
				break;
			case 6: return 173;
				break;
			case 21: return 0;
				break;
			default: return -1;
				break;
		}
	}
};
#endif
#ifndef Jet_h
#define Jet_h
class Jet
{
public:
	Jet(){
		pT=0;
		energy=0;
	}
	Jet(float _pT, float _phi, float _y, float _r){
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
	}
	Jet(double _pT, double _phi, double _y){
		pT=Scalar(_pT);
		phi=Scalar(_phi);
		y=Scalar(_y);
		r=0;
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz){ // calculate eta
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		eta= Scalar(calculateEta(_pT,pz));
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass){ // calculate eta and e
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		energy = Scalar(calculateEnergy(pz));
		eta= Scalar(calculateEta(_pT,pz));
	}
	Jet(float _pT, float _phi, float _y, float _r, float pz, float mass,float energy){ // calculate eta and e
		this->pT =Scalar(_pT);
		this->phi = Scalar(_phi);
		this->y = Scalar(_y);
		this->r = Scalar(_r);
		this->mass = Scalar(mass);
		this->energy = Scalar(energy);
		eta= Scalar(calculateEta(_pT,pz));
	}
	~Jet(){	}
	void setMult(int m){
		mult=m;
	}
	Parton setParton(Parton p1, Parton p2){  
		if (deltaR(p1)<deltaR(p2))//comparision
		{
			parton=p1;	
		}
		else{
			parton=p2;
		}
		return parton;
	}
	inline float deltaPhi(float in){
		float r = TMath::Abs(in-phi.value);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}
	void setParton(Parton p){
		parton=p;
	}
	bool isJetQuark(){
		return parton.getQuark();
	}
	int getmult(){
		return mult;
	}
	Scalar getpT(){
		return pT;
	}
	Scalar getphi(){
		return phi;
	}
	Scalar gety(){
		return y;
	}
	Scalar getr(){
		return r;
	}
	Scalar getEnergy(){
		return energy;
	}
	Scalar operator/(float s){ 
		return pT/s;
	}
	Scalar operator/(Jet s){ 
		return pT/s.getpT();
	}
	bool operator==(Jet s){
		return pT==s.pT&&phi==s.getphi();
	}
	bool operator>(Jet s){
		return pT>s.getpT();
	}
	bool operator<(Jet s){
		return pT<s.getpT();
	}
	bool operator>=(Jet s){
		return pT>=s.getpT();
	}
	bool operator<=(Jet s){
		return pT<=s.getpT();
	}
	bool operator!=(Jet s){
		return pT!=s.getpT();
	}
private:
	Scalar pT=-1;
	Scalar phi =7;
	Scalar y=0;
	Scalar r=-1;
	Scalar eta;
	int mult=0;
	Scalar mass;
	Scalar energy;
	Parton parton;

	float deltaR(Parton p){
	  return TMath::Power((TMath::Power(TMath::Abs(p.geteta()-eta.value),2)+TMath::Power(TMath::Abs(p.getphi()-phi.value),2)),.5);
	}
	float calculateEta(float pt, float pz){
		return .5* TMath::Log((TMath::Power(pt*pt+pz*pz,.5)+pt))/((TMath::Power(pt*pt+pz*pz,.5)-pt));
	}
	float calculateEnergy(float pz){
		return TMath::Power((pT.value)*(pT.value)+pz*pz+mass.value*mass.value,.5);
	}
	
};
#endif
#ifndef DiJet_h
#define DiJet_h
class DiJet
{
public:
	DiJet(Jet j1, Jet j2){
		leading = bigger(j1,j2);
		subleading=smaller(j1,j2);
		calculateR2J2();
	}
	DiJet(Jet j1, Jet j2,bool t){ //for sorted jets 
		if (t)
		{
			isDijet=true;
			leading=j1;
			subleading=j2;
		}
		else{
			leading = bigger(j1,j2);
			subleading=smaller(j1,j2);
		}
		calculateR2J2();
	}
	DiJet(double pt1, double phi1, double pt2,double phi2){
		if (pt1>pt2)
		{
			leading=Jet(pt1,phi1,0,0);
			subleading=Jet(pt2,phi2,0,0);
		}
		else{
			subleading=Jet(pt1,phi1,0,0);
			leading=Jet(pt2,phi2,0,0);
		}
		calculateR2J2();
		makeJetDeltaPhi();
	}
	DiJet(bool f){
		isDijet=f;
	}
	DiJet(){}

	~DiJet(){}
	Jet getleading(){
		return leading;
	}
	Jet getsubleading(){
		return subleading;
	}
	float getR2J2(){
		return r2j2;
	}
	float getDeltaPhi(){
		return jetDeltaPhi;
	}
	void operator=(DiJet d2){
		isDijet=(bool)d2;
		leading=d2.getleading();
		subleading=d2.getsubleading();
		jetDeltaPhi=d2.getDeltaPhi();
	}
	operator bool(){
		return isDijet;
	}
private:
	Jet leading;
	Jet subleading;
	float jetDeltaPhi;
	float photonDeltaPhi;
	float r2j2;
	bool isDijet;

	inline float deltaPhi(float i1, float i2){
		float r = TMath::Abs(i1-i2);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}

	inline float makeJetDeltaPhi(){
		jetDeltaPhi = TMath::Abs(leading.getphi().value-subleading.getphi().value);
		if (jetDeltaPhi>TMath::Pi())
		{
			jetDeltaPhi= 2*TMath::Pi()-jetDeltaPhi;
		}
		return jetDeltaPhi;
	}

	inline void calculateR2J2(){
		r2j2 = (leading.getEnergy().value-subleading.getEnergy().value)/(leading.getEnergy().value+subleading.getEnergy().value);
	}
	
};
#endif
#ifndef myParticle_h
#define myParticle_h
class myParticle
{
public:
	myParticle(){}
	~myParticle(){}
	myParticle(int id, float eT, float phi, float y){
		this->id=id;
		this->eT=eT;
		this->phi=phi;
		this->y=y;
	}
	int getId(){
		return id;
	}
	float geteT(){
		return eT;
	}
	float getphi(){
		return phi;
	}
	float gety(){
		return y;
	}

private:
	int id;
	float eT;
	float phi;
	float y;
};
#endif

inline float deltaPhi(float i1, float i2){
	float r = TMath::Abs(i1-i2);
	if (r>TMath::Pi())
	{
		r= 2*TMath::Pi()-r;
	}
	return r;
}

char genRandomChar()  // Random char generator function.
{
    return alphanum[rand() % 69];
}

string randomString(){
	std::string Str;
	for(unsigned int i = 0; i < 21; ++i)
	{
	    Str += genRandomChar();
	}
	return Str;
}

queue<TH1D*> makeYProjections(TH2F* data, int* binsTodivide, int NHISTS,int* plotcount){ //turns the TH2F into NHISTS TH1F by spliting on the binsTodivide provide upper and lower bin
	queue<TH1D*> r;
	TH1D* temp; 
	for (int i = 0; i < NHISTS; ++i)
	{
		temp = data->ProjectionY(getNextPlotName(plotcount).c_str(),binsTodivide[i],binsTodivide[i+1],"e");
		r.push(temp);
	}
	return r;
}

void drawTemp(TH1* plot){
	TCanvas *tc = new TCanvas();
	plot->Draw();
}

int* getBinsFromValues(TH1* data, bool xaxis,float* values, int SIZE){
	int *bins = new int[SIZE];
	cout<<"Bins:";
	for (int i = 0; i < SIZE; ++i)
	{
		if (xaxis)
		{
			bins[i]=data->GetXaxis()->FindBin(values[i]);
		}
		else{
			bins[i]=data->GetYaxis()->FindBin(values[i]);
		}
		cout<<values[i]<<":"<<bins[i]<<'\n';
	}
	return bins;
}

Scalar getResolution(TH1* plot){
	TF1 *fit = new TF1(randomString().c_str(),"gaus",plot->GetBinContent(1),plot->GetBinContent(plot->GetNbinsX()));
	//plot->Scale(1/plot->Integral());
	plot->Fit(fit);
	float gausData[2];
	gausData[0]= fit->GetParameter(1);
	gausData[1]= fit->GetParameter(2);
	recursiveGaus(plot,fit,gausData,1.5,90);
	Scalar mean(gausData[0],fit->GetParError(1));
	Scalar sigma(gausData[1],fit->GetParError(2));
	//ssanl<<"Mean:"<<gausData[0]<<", "<<gausData[1]<<"="<<gausData[1]/gausData[0]<<'\n';
	return sigma/mean;
}

void plotWithGaus(TH1* plot){
	TCanvas *tc =new TCanvas();
	TF1 *fit = new TF1(randomString().c_str(),"gaus",plot->GetBinContent(1),plot->GetBinContent(plot->GetNbinsX()));
	plot->Scale(1/plot->Integral());
	plot->Fit(fit,"Q0");
	float gausData[2];
	gausData[0]= fit->GetParameter(1);
	gausData[1]= fit->GetParameter(2);
	recursiveGaus(plot,fit,gausData,1.5,90);
	string sigma = "#sigma:"+to_string(gausData[1]);
	plot->Draw();
	fit->SetRange(gausData[0]-5*gausData[1],0,gausData[0]+5*gausData[1],1);
	fit->Draw("same");
	//myText(.2,.3,kBlack,sigma.c_str()); //need to change the build order in root
}

class PlotWithLine
{
public:
	virtual void Draw(){
		main->Draw();
	}
	~PlotWithLine(){ // there might be some mem leakage here 
	}
protected:
	TH1 *main;
	
};

class CutPlot :public PlotWithLine{
public:
	CutPlot(TH1 *main, TLine* cut) : cut(cut){
		this->main  = main;
	}
	~CutPlot(){
		delete cut;
		cut=NULL;
	}
	void Draw(){
		main->Draw();
		cut->Draw("same");
	}
private:
	TLine *cut;
};

class GausPlot :public PlotWithLine
{
public:
	GausPlot(TH1* main, TF1* gaus,double lowBound,double upBound) : gaus(gaus),lowBound(lowBound), upBound(upBound){
		this->main  = main;
	}
	~GausPlot(){
		delete gaus;
		gaus=NULL;
	}
	void Draw(){
		main->Draw();
		gaus->Draw("same");
	}
	double getUpBound(){
		return upBound;
	}

private:
	TF1 *gaus;
	double lowBound;
	double upBound;
};

struct Poly2
{
	Scalar c;
	Scalar b;
	Scalar a;
	Poly2(Scalar _a, Scalar _b, Scalar _c){
		c=_c;
		b=_b;
		a=_a;
	}
	Poly2(){};
};

#ifndef Cluster_h
#define Cluster_h 

class Cluster
{
public:
	Cluster(){}
	Cluster(float pT,float phi, float eta){
		this->pT = pT;
		this->phi = phi;
		this->eta =eta;
	}
	Cluster(float pT,float phi, float eta, int index){
		this->pT = pT;
		this->phi = phi;
		this->eta =eta;
		this->index = index;

	}
	~Cluster(){}
	float setdR(float phi,float eta){
		dR =TMath::Power((TMath::Power(TMath::Abs(eta-this->eta),2)+TMath::Power(deltaPhi(this->phi,phi),2)),.5);
	  	return dR;
	}
	float getdR(){
		return dR;
	}
	float getpT(){
		return pT;
	}
	float geteta(){
		return eta;
	}
	int getIndex(){
		return index;
	}
private:
	float phi;
	float eta;
	float dR;
	float pT;
	int index;
	inline float deltaPhi(float i1, float i2){
		float r = TMath::Abs(i1-i2);
		if (r>TMath::Pi())
		{
			r= 2*TMath::Pi()-r;
		}
		return r;
	}
	
};
#endif

inline float deltaR(float e1, float p1, float e2, float p2){
	  return TMath::Power(TMath::Power(TMath::Abs(e1-e2),2)+TMath::Power(deltaPhi(p1,p2),2),.5);
}

inline float pToE(float x, float y, float z, float mass){
	return quadrature((float)quadrature(x,y),quadrature(z,mass));
}

inline float pToE(TVector3 v, float mass){
	return quadrature((float) quadrature(v.x(),v.y()),(float) quadrature((float)v.z(),mass));
}

/*
class multiTH1F
{
public:
	multiTH1F(){}
	~multiTH1F(){}
	multiTH1F(float max, float min, float plotwidth, int Nbins){
		Nplots = (max-min)/plotwidth;
		for (int i = 0; i < Nplots; ++i)
		{
			v.push_back(new TH1F(getNextPlotName(&plotCount).c_str(),"",Nbins,min+i*plotwidth,min+(i+1)*plotwidth));
			plotmins.push_back(min+i*plotmins);
		}
	}
	void fill(float in){
		v[getPlotN(in)]->Fill(in);
	}
	void normalize(){
		for (std::vector<TH1F*>::iterator i = v.begin(); i != v.end(); ++i)
		{
			(*i)->Scale(1/(*i)->Integral());
		}
	}
	void plot(){
		TCanvas* tc= new TCanvas();
		tc->Divide((int)(v.size()+1)/2,2);
	}

private:
	std::vector<TH1F*> v;
	int Nplots;
	std::vector<float> plotmins;
	int getPlotN(float in){
		if (in<plotmins[0]||in>*plotmins.back())
		{
			return -1;
		}
		int i=1;
		while(in>plotmins[i]){
			i++;
		}
		return i;
	}
	
};*/
