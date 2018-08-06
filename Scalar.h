#ifndef Scalar_h
#define Scalar_h

#include <iostream>
#include <TMath.h>

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
	Scalar operator+(const Scalar &s);
	Scalar operator/(const Scalar &s);
	void operator+=(Scalar s);
	inline void operator*=(Scalar s){
		*this = *this * s;
	}
	inline void operator/=(Scalar s){
		*this = *this / s;
	}
	inline void operator-=(Scalar s){
		*this = *this - s;
	}
	inline Scalar operator+(float s){
		return Scalar(this->value+s);
	}
	Scalar operator-(Scalar s);
	inline Scalar operator-(float s){
		return Scalar(this->value-s);
	}
	inline void operator=(Scalar s){
		value=s.value;
		uncertainty=s.uncertainty;
	}
	inline bool operator==(Scalar s){
		return value==s.value;
	}
	inline bool operator==(float s){
		return value==s;
	}
	inline bool operator!=(float s){
		return value!=s;
	}
	inline bool operator!=(double s){
		return value!= (float) s;
	}
	inline bool operator>(Scalar s){
		return value>s.value;
	}
	inline bool operator<(Scalar s){
		return value<s.value;
	}
	inline bool operator>=(Scalar s){
		return value>=s.value;
	}
	inline bool operator<=(Scalar s){
		return value<=s.value;
	}
	inline bool operator!=(Scalar s){
		return value!=s.value;
	}
	Scalar operator*(Scalar s);
	Scalar operator*(float s);
	Scalar operator*(double s);
	Scalar operator/(float s);
	Scalar pow(double n);
	Scalar log(float base);

	/*Scalar totalaverage(queue<Scalar> ins){
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
	}*/
	inline operator float(){return value;}
	inline operator double(){return (double) value;}
	inline friend std::ostream& operator<<(std::ostream& os, Scalar const & tc) {
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