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
