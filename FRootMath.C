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
