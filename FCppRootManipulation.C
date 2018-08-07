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
	queue<T> fillQueue(T in, const int SIZE){
		queue<T> out;
		for (int i = 0; i < SIZE; ++i)
		{
			out.push(in);
		}
		return out;
	}
