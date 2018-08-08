#include "Scalar.h"
#include "TMath.h"
#include <iostream>
using namespace std;

//not sure of all the inheritence rules on virtualization/abstracttion want it to do all the same things as a Scalar but check for <PI b4 all the returns
class Angle : public Scalar //all the members of Scalar become members of Angle with equal scope
{
public:
	Angle(){
		value=uncertainty=0;
	}
	~Angle(){}
	Angle(float value){
		this->value=value;
		uncertainty=0;
	}
	Angle(float value, float uncertainty){
		this->value=value;
		this->uncertainty=uncertainty;
	}
	Angle(const Scalar &s){
		this->value=s.value;
		uncertainty=s.uncertainty;
	}
	Angle operator+(float s){
		Angle next;
		next.value = this->value + s;
		return next;
	}
	Angle operator-(Scalar s){
		Angle next;
		next.value = this->value - s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		return next;
	}
	Angle operator-(float s){
		Angle next;
		next.value = this->value - s;
		return next;
	}
	Angle operator*(Scalar s){
		Angle next;
		next.uncertainty = (value*s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		next.value = value*s.value;
		makeIn2Pi();
		return next;
	}
	Angle operator*(float s){
		Angle next;
		next.value = value*s;
		next.uncertainty = (value*s)*uncertainty/value;
		if (value==0||s==0)
		{
			next.uncertainty=0;
		}
		makeIn2Pi();
		return next;
	}
	Angle operator/(float s){ //progate uncertainty
		Angle next;
		next.uncertainty=(value/s)*uncertainty/value;
		if (value==0)
		{
			next.uncertainty=0;
		}
		next.value = value/s;
		makeIn2Pi();
		return next;
	}
	Angle operator+(const Scalar &s){
		//cout<<"Addition "<<value<<" + "<<s.value<<"\n";;
		Angle next = Scalar();
		next.value = this->value + s.value;
		next.uncertainty = quadrature(this->uncertainty, s.uncertainty);
		makeIn2Pi();
		return next;
	}
	Angle operator/(const Scalar &s){
		Angle next = Scalar();
		next.uncertainty = (value/s.value)*quadrature(this->uncertainty/value, s.uncertainty/s.value);
		next.value = value/s.value;
		if (value==0||s.value==0)
		{
			next.uncertainty=0;
		}
		makeIn2Pi();
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
	void operator=(Scalar s){
		value=s.value;
		uncertainty=s.uncertainty;
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
	Angle pow(double n){
		Angle r;
		if (n==0)
		{
			r.value=1;
			r.uncertainty=0;
			return r;
		}
		r.value = TMath::Power((double)value,n);
		r.uncertainty=(r.value*n*uncertainty/value);
		makeIn2Pi();
		return r;
	}
	Angle log(float base){
		Angle r;
		r.value = TMath::Log(value)/TMath::Log(base);
		r.uncertainty = uncertainty/(value*TMath::Log(base));
		makeIn2Pi();
		return r;
	}
	friend std::ostream& operator<<(std::ostream& os, Angle const & tc) {
       return os <<"Radian:" << tc.value <<char(241)<<tc.uncertainty<<'\n';
    }
private:
	inline void makeIn2Pi(){
		if(value<0) value+=2*TMath::Pi();
	}
};

template<class T>
T quadrature(T d1, T d2){
  return sqrt((double)d1*d1+d2*d2);
}
