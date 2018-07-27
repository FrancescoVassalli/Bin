#include "TLegend.h"
#include "TH1F.h"
#include <limits.h>
#include <queue>
#include <iostream>
class Scalar;

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
	template<class T>
	void arrayMultiply(T *a, T* m,int SIZE){
		for (int i = 0; i < SIZE; ++i)
		{
			a[i] = a[i]*m[i];
		}
	}
	template<class T>
	T sum(T *a, int SIZE){
		T sum=0;
		for (int i = 0; i < SIZE; ++i)
		{
			sum+=a[i];
		}
	}
	template<class T>
	T sum(queue<T> q){
		//queue<T> q = clone(a);
		T sum=0;
		while(!q.empty()){
			sum+=q.front();
			q.pop();
		}
		return sum;
	}
	template<class T>
	T quadrature(T d1, T d2){
		return TMath::Sqrt((double)d1*d1+d2*d2);
	}
	template<class T>
	T quadrature(T* a, int SIZE){
		T* b = clone(a);
		arrayMultiply(b,b,SIZE);
		T r = sum(b,SIZE);
		return TMath::Sqrt(r);

	}
	template<class T>
	void divideArray(int n, T* a, T d){
		for (int i = 0; i < n; ++i)
		{
			a[i] /=d;
			//cout<<a[i]<<'\n';
		}
	}
	template<class T>
	bool inRange(T in, T low, T high){
		return in>=low && in<=high;
	}
	template<class T>
	void divideArray(int n, T* a, T* d){
		for (int i = 0; i < n; ++i)
		{
			a[i] /=d[i];

		}
	}
	template<class T>
	T errorDivide(T d1, T e1, T d2, T e2){
		T r = d1/d2;
		return r * quadrature(e1/d1,e2/d2);
	}
	template<class T>
	T errorDivide(T d1, T e1, int d2){
		T r = d1/d2;
		return r * e1/d1;
	}
	template<class T>
	void errorDivideArray(int n,T* values, T* errors, T divisor, T divisorError){
		for (int i = 0; i < n; ++i)
		{
			errors[i] = errorDivide(values[i],errors[i],divisor,divisorError);
		}
		divideArray(n,values,divisor);
	}
	template<class T>
	void errorDivideArray(int n,T* values, T* errors, T* divisor, T* divisorError){
		for (int i = 0; i < n; ++i)
		{
			errors[i] = errorDivide(values[i],errors[i],divisor[i],divisorError[i]);
		}
		divideArray(n,values,divisor);
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
	template<class T>
	T* queueToArray(queue<T> q){
		T* r = new T[q.size()];
		int i=0;
		while(!q.empty())
		{
			r[i++] = q.front();
			q.pop();
		}
		return r;
	}
	/*returns every index i in an array where the i value is unequal to the value at i-1 always starts with 0 and ends with SIZE*/
	template<class T>
	vector<int>* sameValueIndices(const int SIZE, T* a){
		vector<int> *maxes = new vector<int>(0);
		maxes->push_back(0);
		if(SIZE<=1){
			maxes->push_back(0);
			return maxes;
		}
		T temp;
		temp = a[0];
		for (int i = 1; i < SIZE; ++i)
		{
			if(a[i]!=temp){
				maxes->push_back(i);
				temp=a[i];
				//cout<<i<<endl;
			}
		}
		maxes->push_back(SIZE);
		return maxes;
	}
	/* pass a vector of the indices from an array you want*/
	template<class T>
	T* valuesAt(T* a,std::vector<int> v){
		T* r = new T[v.size()];
		for (int i = 0; i < v.size(); ++i)
		{
			r[i] = a[v[i]];
		}
		return r;
	}
	/*inclusive*/
	template<class T>
	T* partialArray(T* a, int start, int end){
		T *r = new T[end-start+1];
		int count=0;
		for (int i = start; i <= end; ++i)
		{
			r[count++]=a[i];
		}
		return r;
	}
	template<class T>
	T* combineArray(T* a, int s1, T* b,int s2){
		T* r = new T[s1+s1];
		for (int i = 0; i < s1; ++i)
		{
			r[i]=a[i];
		}
		for (int i = 0; i < s2; ++i)
		{
			r[i+s1] =b[i];
		}
		return r;
	}
	/*exludes index a and b*/
	template<class T>
	T* removeIndex(T* in, int a, int b, int *SIZE){
		T *r1 = partialArray(in,0,a-1);
		T *r2 = partialArray(in,b+1,*SIZE-1);
		T *rf = combineArray(r1,a,r2,*SIZE-b);
		*SIZE = *SIZE-(b-a)-1;
		return rf;
	}
	template<class T>
	queue<T> vectorToQ(std::vector<T> v){
		queue<T> r;
		for (typename std::vector<T>::iterator i = v.begin(); i != v.end(); ++i)
		{
			r.push(*i);
		}
		return r;
	}
	/*exludes index a and b*/
	template<class T>
	queue<T> removeIndex(queue<T> in, int a, int b){
		std::vector<T> v = queueToVector(in);
		int count =0;
		for (typename std::vector<T>::iterator i = v.begin(); i != v.end(); ++i)
		{
			if (count>=a&&count<=b)
			{
				v.erase(i);
			}
			count++;
		}
		queue<T> r = vectorToQ(v);
		return r;
	}
	/*template<class T>
	void removeIndex(queue<queue<T>> *in, int a, int b){
		queue<queue<T>> *r = new queue<queue<T>>;
		int count=0;
		while(!in->empty()){
			if (count<a&&count>b)
			{
				r->push(in->front());			
			}
			in->pop();
			count++;
		}
		delete in;
		in=r;
	}*/
	template<class T>
	std::vector<T> queueToVector(queue<T> q){
		std::vector<T> v;
		while(!q.empty()){
			v.push_back(q.front());
			q.pop();
		}
		return v;
	}
	template<class T>
	T bigger(T x, T y){
		if (x>y)
		{
			return x;
		}
		else{
			return y;
		}
	}
	template<class T>
	T smaller(T x, T y){
		if (x<y)
		{
			return x;
		}
		else{
			return y;
		}
	}
	/* queues a-b inclusive are put into queue a and the other queues are removed */
	template<class T>
	queue<queue<T>> mergeQueues(queue<queue<T>> in, int a, int b){
		std::vector<queue<T>> v = queueToVector(in);
		queue<queue<T>> out;
		int count=0;
		queue<T> place = v.at(a);
		for (typename std::vector<queue<T>>::iterator i = v.begin(); i != v.end(); ++i)
		{
			if (count<a||count>b)
			{
				//cout<<"normal push"<<'\n';
				out.push(*i);
			}
			else{
				if (count!=a)
				{
					//cout<<"place push"<<'\n';
					while(!i->empty()){
							place.push(i->front());
							i->pop();
						}
					if (count==b)
					{
						//cout<<"final push"<<'\n';
						out.push(place);
					}
				}
			}
			count++;
		}
		return out;
	}
	template<class T>
	queue<int> arrayNonZero(T* a, int SIZE){
		queue<int> r;
		for (int i = 0; i < SIZE; ++i)
		{
			if (a[i]!=0)
			{
				r.push(i);
			}
		}
		return r;
	}
	/*inclusive*/
	template<class T>
	queue<T> arrayToQueue(T* a, int start, int end){
		queue<T> r;
		while(start<=end){
			r.push(a[start]);
		}
		return r;
	}
	template<class T>
	queue<T> arrayToQueue(const int SIZE, T* a){
		queue<T> r;
		for(int i=0; i<SIZE;++i){
			r.push(a[i]);
		}
		return r;
	}
	/** returns a 2-D array where each subarray is the parts of a broken down by the output of maxes and partial Array*/
	/*template<class T>
	queue<queue<T>> maxBrokenArray(const int SIZE, T* a){
		vector<int> *maxis = maxes(SIZE,a);
		for (unsigned int i = 0; i < maxis->size(); ++i)
		{
			cout<<maxis->at(i)<<'\n';
		}
		cout<<endl;
		queue<queue<T>> out;
		for (int i = 0; i < maxis->size()-1; ++i)
		{
			out.push(arrayToQueue(a,(*maxis)[i],(*maxis)[i+1]));
		}
		return out;
	}*/
	template<class T>
	T* sqrtArray(int n, T* in){
		T* t= new T[n];
		for (int i = 0; i < n; ++i)
		{
			t[i] = TMath::Sqrt(in[i]);
		}
		return t;
	}
	template<class T>
	T* sigmaNtoUncertainty(int SIZE, T* sigma, T* n, T* sigmaerror){
		const float factor = 1.465; // the inverse of 0.6827 which is the percentage of the area under a gaussian when integrated over +- 1 sigma
		T *r = new T[SIZE];
		T *r2 = new T[SIZE];
		for (int i = 0; i < SIZE; ++i)
		{
			r[i]= factor*sigma[i];
			r2[i] = sigmaerror[i];
		}
		errorDivideArray(SIZE,r,r2,sqrtArray(SIZE,n),zeroArray(SIZE,n));
		return r;
	}
	template<class T>
	T max(int SIZE, T* x){
		T temp=0;
		for (int i = 0; i < SIZE; ++i)
		{
			if (x[i]>temp)
			{
				temp=x[i];
			}
		}
		return temp;
	}
	template<class T>
	T max(queue<T> in){
		T temp=0;
		while(!in.empty()){
			if(in.front()>=temp){
				temp=in.front();
			}
			in.pop();
		}
		return temp;
	}
	template<class T>
	T min(queue<T> in){
		T temp=-1;
		while(!in.empty()){
			if(in.front()<=temp||temp==-1){
				temp=in.front();
			}
			in.pop();
		}
		return temp;
	}
	/*template<class T>
	T min(T t1, T t2){
		if (t1<t2)
		{
			return t1;
		}
		else return t2;
	}*/
	template<class T>
	T average(int SIZE, T* x){
		T sum=0;
		for (int i = 0; i < SIZE; ++i)
		{
			sum= sum +x[i];
		}
		return sum/SIZE;
	}
	template<class T>
	T average(queue<T> in){
		T sum=0;
		const int SIZE = in.size();
		while(!in.empty()){
			sum= sum+in.front();
			in.pop();
		}
		//cout<<sum<<"\n";
		return sum/SIZE;
	}
	/*float systematicError(const int SIZE, float* means){ // by extreme - mean over sqrt(3)
		return (max<float>(SIZE,means)-average<float>(SIZE,means))/TMath::Sqrt(3);
	}*/
	template<class T>
	T systematicError(queue<T> means){ // by extreme 
		if (means.size()<2)
		{
			return 0;
		}
		T top = max<T>(means);
		T bottom = min<T>(means);
		T a = average<T>(means);
		if (top-a>bottom-a)
		{
			return((top-a)/TMath::Sqrt(3));
		}
		else{
			return ((bottom-a)/TMath::Sqrt(3));
		}
	}
	template<class T>
	queue<T> groupSystematic(queue<queue<T>> groups){
		queue<T> out;
		while(!groups.empty()){
			out.push(systematicError(groups.front()));
			groups.pop();
		}
		return out;
	}
	template<class T>
	void oleSwitcheroo(T* xp, T* yp)
	{
    	T temp = *xp;
    	*xp = *yp;
    	*yp = temp;
	}
	template<class T>
	queue<queue<T>> breakArray(T* a, std::vector<int> breaks){
		queue<queue<T>> r;
		int breakcounter=0;
		while(breaks.size()-breakcounter>1){
			//cout<<breaks.at(breakcounter)<<" and "<<breaks.at(breakcounter+1)<<endl;
			queue<T> temp;
			for (int i = breaks.at(breakcounter); i < breaks.at(breakcounter+1); ++i)
			{
				//cout<<a[i]<<endl;
				temp.push(a[i]);
			}
			r.push(temp);
			breakcounter++;
			//cout<<r.back().front()<<'\n';
		}
		return r;
	}
	template<class T>
	queue<T> averageList(queue<queue<T>> q){
		queue<T> r;
		while(!q.empty()){
			r.push(average(q.front()));
			q.pop();
		}
		return r;

	}
	void recursiveGaus(TH1* h, TF1* gaus, float* data, float sigmadistance,int lazyMan=0){
	    h->Fit(gaus,"Q0","",data[0]-sigmadistance*data[1],data[0]+sigmadistance*data[1]);
	    if(data[0]!=gaus->GetParameter(1)){
	    	if(lazyMan == 100){return;}
	        data[0] = gaus->GetParameter(1);
	        data[1] = gaus->GetParameter(2);
	        lazyMan++;
	        recursiveGaus(h,gaus,data,sigmadistance,lazyMan);
	    }
	    else{
	        data[0] = gaus->GetParameter(1);
	        data[1] = gaus->GetParameter(2);
	        return;
	    }
	}
	/* note that SIZE will change to the size of the returned array*/
	template<class T>
	T* yAveragex(unsigned int* SIZE, T* y, T* x){
		queue<queue<T>> broken = breakArray(y,*sameValueIndices(*SIZE,x));
		queue<T> averages;
		while(!broken.empty()){
			averages.push(average(broken.front()));
			broken.pop();
		}
		T* r= queueToArray(averages);
		*SIZE = averages.size();
		return r;
	}
	template<class T>
	void averageYXPropogate(unsigned int* SIZE, T* y, T* x, T*dy){
		queue<queue<T>> broken = breakArray(y,*sameValueIndices(*SIZE,x));
		queue<queue<T>> brokenUncertain = breakArray(dy,*sameValueIndices(*SIZE,x));
		queue<T> averages;
		queue<T> uncertainty;
		while(!broken.empty()){
			int count  = broken.front().size();
			T temp = sum(broken.front());
			T uTemp = sum(brokenUncertain.front());
			uTemp=errorDivide(temp,uTemp,count);
			temp/=count;
			averages.push(temp);
			uncertainty.push(uTemp);
			broken.pop();
			brokenUncertain.pop();
		}
		y= queueToArray(averages);
		dy= queueToArray(uncertainty);
		*SIZE = averages.size();
	}

	template<class T>
	T* uniqueSortedArrayValues(int SIZE, T* a){
		std::vector<int> *diffs = sameValueIndices(SIZE,a);
		T* r = new T[diffs->size()];
		for(unsigned int i=0; i<diffs->size();i++){
			r[i] = a[diffs->at(i)];
		}
		return r;
	}
	template<class T>
	queue<T> fillQueue(T in, const int SIZE){
		queue<T> out;
		for (int i = 0; i < SIZE; ++i)
		{
			out.push(in);
		}
		return out;
	}
	
#ifndef Scalar_h
#define Scalar_h
class Scalar
{
public:
	Scalar(){
		value=uncertainty=0;
	}
	Scalar(float value){
		this->value=value;
		uncertainty=0;
	}
	Scalar(float value, float uncertainty){
		this->value=value;
		this->uncertainty=uncertainty;
	}
	Scalar(const Scalar &s){
		this->value=s.value;
		uncertainty=s.uncertainty;
	}
	~Scalar(){}; 
	//operators might need to return class types but I'm not sure
	Scalar operator+(const Scalar &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		Scalar next = Scalar();
		next.value = this->value + s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Scalar operator/(const Scalar &s){
		Scalar next = Scalar();
		next.uncertainty = (value/s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		next.value = value/s.value;
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		return next;
	}
	void operator+=(Scalar s){
		value+=s.value;
		uncertainty = quadrature(this->uncertainty, s.uncertainty);
	}
	void operator*=(Scalar s){
		*this = *this * s;
	}
	void operator/=(Scalar s){
		*this = *this / s;
	}
	void operator-=(Scalar s){
		*this = *this - s;
	}
	Scalar operator+(float s){
		Scalar next;
		next.value = this->value + s;
		return next;
	}
	Scalar operator-(Scalar s){
		Scalar next;
		next.value = this->value - s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Scalar operator-(float s){
		Scalar next;
		next.value = this->value - s;
		return next;
	}
	void operator=(Scalar s){
		value=s.value;
		uncertainty=s.uncertainty;
	}
	bool operator==(Scalar s){
		return value==s.value;
	}
	bool operator==(float s){
		return value==s;
	}
	bool operator!=(float s){
		return value!=s;
	}
	bool operator!=(double s){
		return value!= (float) s;
	}
	bool operator>(Scalar s){
		return value>s.value;
	}
	bool operator<(Scalar s){
		return value<s.value;
	}
	bool operator>=(Scalar s){
		return value>=s.value;
	}
	bool operator<=(Scalar s){
		return value<=s.value;
	}
	bool operator!=(Scalar s){
		return value!=s.value;
	}
	Scalar operator*(Scalar s){
		Scalar next;
		next.uncertainty = (value*s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		next.value = value*s.value;
		return next;
	}
	Scalar operator*(float s){
		Scalar next;
		next.value = value*s;
		next.uncertainty = (value*s)*uncertainty/value;
		if (value==0||s==0)
		{
			next.uncertainty=0;
		}
		return next;
	}
	Scalar operator*(double s){
		Scalar next;
		next.value = value*(float)s;
		next.uncertainty = (value*(float)s)*uncertainty/value;
		if (value==0||s==0)
		{
			next.uncertainty=0;
		}
		return next;
	}
	Scalar operator/(float s){ //progate uncertainty
		Scalar next;
		next.uncertainty=(value/s)*uncertainty/value;
		if (value==0)
		{
			next.uncertainty=0;
		}
		next.value = value/s;
		return next;
	}
	Scalar pow(double n){
		Scalar r;
		if (n==0)
		{
			r.value=1;
			r.uncertainty=0;
			return r;
		}
		r.value = TMath::Power((double)value,n);
		r.uncertainty=(r.value*n*uncertainty/value);
		return r;
	}
	Scalar log(float base){
		Scalar r;
		r.value = TMath::Log(value)/TMath::Log(base);
		r.uncertainty = uncertainty/(value*TMath::Log(base));
		return r;
	}

	Scalar totalaverage(queue<Scalar> ins){
		Scalar temp=0;
		float numerator=0;
		float denominator=0;
		while(!ins.empty()){
			temp+=ins.front();
			numerator+=ins.front().value/(ins.front().uncertainty*ins.front().uncertainty);
			denominator+=1/(ins.front().uncertainty*ins.front().uncertainty);
			ins.pop();
		}
		temp.value=numerator/denominator;
		return temp;
	}
	operator float(){return value;}
	operator double(){return (double) value;}
	friend std::ostream& operator<<(std::ostream& os, Scalar const & tc) {
       return os <<"Scalar:" << tc.value <<char(241)<<tc.uncertainty<<'\n';
    }

	float value;
	float uncertainty;
protected:
	template<class T>
T quadrature(T d1, T d2){
	return TMath::Sqrt((double)d1*d1+d2*d2);
}

template<class T>
T quadrature(T* a, int SIZE){
	T* b = clone(a);
	arrayMultiply(b,b,SIZE);
	T r = sum(b,SIZE);
	return TMath::Sqrt(r);
}
	
};
#endif
#ifndef ScalarAsymmetric_h
#define ScalarAsymmetric_h
class ScalarAsymmetric : public Scalar
{
public:
	ScalarAsymmetric(){
		value=uncertaintyUp=uncertaintyDown=0;
	}
	ScalarAsymmetric(float value){
		this->value=value;
		uncertaintyUp=uncertaintyDown=0;
	}
	ScalarAsymmetric(float value, float uncertainty){
		this->value=value;
		uncertaintyUp=uncertaintyDown=uncertainty;
	}
	ScalarAsymmetric(float value, float uncertaintyUp,float uncertaintyDown){
		this->value=value;
		this->uncertaintyUp=uncertaintyUp;
		this->uncertaintyDown=uncertaintyDown;
	}
	~ScalarAsymmetric(){}; 

	ScalarAsymmetric operator+(const ScalarAsymmetric &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		ScalarAsymmetric next;
		next.value = this->value + s.value;
		next.uncertaintyUp = quadrature(this->uncertaintyUp, s.uncertaintyUp);
		next.uncertaintyDown = quadrature(this->uncertaintyDown, s.uncertaintyDown);
		return next;
	}
	ScalarAsymmetric operator+(const Scalar &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		ScalarAsymmetric next;
		next.value = this->value + s.value;
		next.uncertaintyUp = quadrature(this->uncertaintyUp, s.uncertainty);
		next.uncertaintyDown = quadrature(this->uncertaintyDown, s.uncertainty);
		return next;
	}
	ScalarAsymmetric operator/(const Scalar &s){
		ScalarAsymmetric next;
		next.uncertaintyUp = (value/s.value)*quadrature(this->uncertaintyUp/value, s.uncertainty/s.value);
		next.uncertaintyDown = (value/s.value)*quadrature(this->uncertaintyDown/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		next.value = value/s.value;
		return next;
	}
	ScalarAsymmetric operator/(const ScalarAsymmetric &s){
		ScalarAsymmetric next;
		next.uncertaintyUp = (value/s.value)*quadrature(this->uncertaintyUp/value, s.uncertaintyUp/s.value);
		next.uncertaintyDown = (value/s.value)*quadrature(this->uncertaintyDown/value, s.uncertaintyDown/s.value);
		if (value==0||s.value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		next.value = value/s.value;
		return next;
	}
	void operator+=(Scalar s){
		*this = *this+s;
	}
	void operator+=(ScalarAsymmetric s){
		*this = *this+s;
	}
	void operator*=(Scalar s){
		*this = *this * s;
	}
	void operator/=(Scalar s){
		*this = *this / s;
	}
	void operator-=(Scalar s){
		*this = *this - s;
	}
	void operator*=(ScalarAsymmetric s){
		*this = *this * s;
	}
	void operator/=(ScalarAsymmetric s){
		*this = *this / s;
	}
	void operator-=(ScalarAsymmetric s){
		*this = *this - s;
	}
	ScalarAsymmetric operator+(float s){
		ScalarAsymmetric next;
		next.value = this->value + s;
		return next;
	}
	ScalarAsymmetric operator-(Scalar s){
		ScalarAsymmetric next;
		next.value = this->value - s.value;
		next.uncertaintyUp = quadrature(this->uncertaintyUp, s.uncertainty);
		next.uncertaintyDown = quadrature(this->uncertaintyDown, s.uncertainty);
		return next;
	}
	ScalarAsymmetric operator-(ScalarAsymmetric s){
		ScalarAsymmetric next;
		next.value = this->value - s.value;
		next.uncertaintyUp = quadrature(this->uncertaintyUp, s.uncertaintyUp);
		next.uncertaintyDown = quadrature(this->uncertaintyDown, s.uncertaintyDown);
		return next;
	}
	ScalarAsymmetric operator-(float s){
		ScalarAsymmetric next;
		next.value = this->value - s;
		return next;
	}
	void operator=(Scalar s){
		value=s.value;
		uncertaintyUp=uncertaintyDown=s.uncertainty;
	}
	void operator=(ScalarAsymmetric s){
		value=s.value;
		uncertaintyUp=s.uncertaintyUp;
		uncertaintyDown=s.uncertaintyDown;
	}
	bool operator==(Scalar s){
		return value==s.value;
	}
	bool operator>(Scalar s){
		return value>s.value;
	}
	bool operator<(Scalar s){
		return value<s.value;
	}
	bool operator>=(Scalar s){
		return value>=s.value;
	}
	bool operator<=(Scalar s){
		return value<=s.value;
	}
	bool operator!=(Scalar s){
		return value!=s.value;
	}
	ScalarAsymmetric operator*(ScalarAsymmetric s){
		ScalarAsymmetric next;
		next.uncertaintyUp = (value*s.value)*quadrature(this->uncertaintyUp/value, s.uncertaintyUp/s.value);
		next.uncertaintyDown = (value*s.value)*quadrature(this->uncertaintyDown/value, s.uncertaintyDown/s.value);
		next.value = value*s.value;
		if (value==0||s.value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		return next;
	}
	ScalarAsymmetric operator*(Scalar s){
		ScalarAsymmetric next;
		next.uncertaintyUp = next.uncertaintyDown= (value*s.value)*quadrature(this->uncertaintyUp/value, s.uncertainty/s.value);
		next.uncertaintyDown = (value*s.value)*quadrature(this->uncertaintyDown/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		next.value = value*s.value;
		return next;
	}
	ScalarAsymmetric operator*(float s){
		ScalarAsymmetric next;
		next.value = value*s;
		next.uncertaintyUp = (value*s)*uncertaintyUp/value;
		next.uncertaintyDown = (value*s)*uncertaintyDown/value;
		if (value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		return next;
	}
	ScalarAsymmetric operator/(float s){ //progate uncertainty
		ScalarAsymmetric next;
		next.uncertaintyUp=(value/s)*uncertaintyUp/value;
		next.uncertaintyDown=(value/s)*uncertaintyDown/value;
		next.value = value/s;
		if (value==0)
		{
			next.uncertaintyUp=next.uncertaintyDown=0;
		}
		return next;
	}
	ScalarAsymmetric pow(double n){
		ScalarAsymmetric r;
		if (n==0)
		{
			r.value=1;
			r.uncertainty=0;
			return r;
		}
		r.value = TMath::Power((double)value,n);
		r.uncertaintyUp=(r.value*n*uncertaintyUp/value);
		r.uncertaintyDown=(r.value*n*uncertaintyDown/value);
		return r;
	}
	ScalarAsymmetric log(float base){
		ScalarAsymmetric r;
		r.value = TMath::Log(value)/TMath::Log(base);
		r.uncertaintyUp = uncertaintyUp/(value*TMath::Log(base));
		r.uncertaintyDown = uncertaintyDown/(value*TMath::Log(base));
		return r;
	}
	void setSymmetric(float uncertainty){
		uncertaintyUp=uncertaintyDown=uncertainty;
	}
	friend ostream& operator<<(ostream& os, ScalarAsymmetric const & tc) {
        return os <<"Scalar: " << tc.value <<"+"<<tc.uncertaintyUp<<"-"<<tc.uncertaintyDown<<'\n';
    }

float value;
float uncertaintyUp;
float uncertaintyDown;
	
};
#endif
#ifndef Point_h
#define Point_h
struct Point
{
	Scalar x;
	Scalar y;
	Point (Scalar _x, Scalar _y){
		x=_x;
		y=_y;
	}
	Point(){}
};
#endif

template<class T>
struct Pair
{
	T x;
	T y;
	Pair();
	Pair(T _x, T _y){
		x=_x;
		y=_y;
	}
};

//probably memory leaks 
queue<queue<Point>> groupPointsByX(Point* a,int *SIZE){
	queue<queue<Point>> r;
	queue<float> interest;
	/*cout<<"check"<<'\n';
	for (int i = 0; i < *SIZE; ++i)
	{
		cout<<a[i].y.value<<'\n';
	}*/
	int groupCount=0;
	//cout<<"Grouping:"<<'\n';
	Scalar temp = a[0].x;
	queue<Point> tempQ;
	for (int i = 0; i < *SIZE; ++i)
	{	
		//cout<<"Group"<<groupCount<<'\n';
		if (a[i].x==temp)
		{
			//cout<<a[i].x.value<<": "<<a[i].y.value<<'\n';
			tempQ.push(a[i]);
		}
		else{
			r.push(tempQ);
			tempQ=queue<Point>();
			temp=a[i].x;
			i-=1;
			groupCount++;
		}
	}
	r.push(tempQ);
	*SIZE = r.size();
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
#ifndef Photon_h
#define Photon_h
class Photon
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
	Photon(double _pT,double _phi, double _eta, bool process, std::queue<myParticle> all){
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
	}
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
	void setParton(Parton p){
		parton=p;
	}
	bool isDirect()const{
		return direct;
	}
	int getPosition()const{
		return position;
	}
	float findIsoEt(queue<myParticle> all){
		isoEt=0;
		while(!all.empty()){
			if (inCone(all.front().gety(),all.front().getphi()))
			{
				isoEt+=all.front().geteT();
			}
			all.pop();
		}
		return isoEt;
	}

	float getIsoEt()const{
		return isoEt;
	}
	
	Point getAngle()const{
		Point r;
		r.x=phi;
		r.y=eta;
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
	Parton parton;
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