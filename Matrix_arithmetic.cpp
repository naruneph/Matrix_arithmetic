// Variant 17 

#include <iostream>
#include <cmath>  
#include <string>  

//Very useful comment

using namespace std;


class RMatrix{
	double **p;
	int row,col;

	void creator();
	double determ(double**,int);
	void cofactor(double **,int ,double **);

public:
	RMatrix();
	RMatrix(int,int);
	RMatrix(const RMatrix &);
	~RMatrix();

	RMatrix & operator =(const RMatrix &);

	friend	ostream & operator <<(ostream &, const RMatrix&);
	friend 	istream & operator >>(istream &, RMatrix &);

	RMatrix & operator +=(const RMatrix &);
	RMatrix & operator -=(const RMatrix &);
	RMatrix & operator *=(const RMatrix &);
	RMatrix & operator *=(double a);
	RMatrix & operator /=(double a);

	bool operator ==(const RMatrix &)const;
	bool operator !=(const RMatrix &)const;

	bool is_square()const;
	bool is_diagonal()const;
	bool is_zero()const;
	bool is_unit()const;
	bool is_symmetric()const;
	bool is_upper_triangular()const;
	bool is_lower_triangular()const;

	RMatrix transpose();
	double determinant();
	RMatrix inverse();
	RMatrix power(int n);

	double norm()const;

};



void RMatrix::creator(){
	p = new double*[row];
	for (int i=0;i<row; ++i){
		p[i] = new double[col];
	}
}

RMatrix::RMatrix(): row(1),col(1){
	creator();
	p[0][0]=0;

}

RMatrix::RMatrix(int x,int y): row(x),col(y){
	creator();
	for (int i=0; i<row; ++i){
		for (int j=0; j<col; ++j){
			p[i][j]=0;
		}
	}
}

RMatrix::RMatrix(const RMatrix & a): row(a.row), col(a.col){
	creator();
	for (int i=0; i<row; ++i){
		for (int j=0; j<col; ++j){
			p[i][j]=a.p[i][j];
		}
	}

}

RMatrix::~RMatrix(){
	for (int i =0; i<row; ++i){
		delete[] p[i];
	}
	delete[] p;
}


RMatrix & RMatrix::operator=(const RMatrix & a){
	if (this == &a){
		return *this;
	}

	if (row != a.row || col != a.col){
		for (int i =0; i<row; ++i){
			delete[] p[i];
		}
		delete[] p;

		row=a.row;
		col=a.col;
		creator();
	}

	for (int i=0; i<row; ++i){
		for (int j=0; j<col; ++j){
			p[i][j]=a.p[i][j];
		}
	}
	return *this;
}

ostream & operator <<(ostream & out, const RMatrix & a){
	for (int i = 0; i<a.row; ++i){
		out << a.p[i][0];
		for (int j = 1; j<a.col; ++j){
			out << " " << a.p[i][j];
		}
		out << endl;
	}
	return out;
}

istream & operator >>(istream & in, RMatrix & a){
	for (int i = 0; i < a.row; ++i){
        for (int j = 0; j < a.col; ++j) {
            in >> a.p[i][j];
        }
    }
    return in;

}

RMatrix & RMatrix::operator +=(const RMatrix & a){
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            p[i][j] += a.p[i][j];
        }
    }
    return *this;
}

RMatrix & RMatrix::operator -=(const RMatrix & a){
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            p[i][j] -= a.p[i][j];
        }
    }
    return *this;
}

RMatrix & RMatrix::operator *=(const RMatrix & a){
	if(this->col == a.row){

		RMatrix temp(row, a.col);

		for (int i = 0; i < temp.row; ++i) {
        	for (int j = 0; j < temp.col; ++j) {
            	for (int k = 0; k < col; ++k) {
                temp.p[i][j] += (p[i][k] * a.p[k][j]);
            	}
        	}
        }
        *this = temp;
        return *this;
	}
//	else{}
}

RMatrix & RMatrix::operator *=(double a){
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            p[i][j] *= a;
        }
    }
    return *this;
}

RMatrix & RMatrix::operator /=(double a){
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            p[i][j] /= a;
        }
    }
    return *this;
}


RMatrix operator +(const RMatrix & a,const RMatrix & b){
	RMatrix temp(a);
	return(temp += b);
}

RMatrix operator -(const RMatrix & a,const RMatrix & b){
	RMatrix temp(a);
	return(temp -= b);
	
}
RMatrix operator *(const RMatrix & a,const RMatrix & b){
	RMatrix temp(a);
	return(temp *= b);	
}
RMatrix operator *(const RMatrix & a,double  b){
	RMatrix temp(a);
	return(temp *= b);	
}
RMatrix operator *(double  b,const RMatrix & a){
	RMatrix temp(a);
	return(temp *= b);
}
RMatrix operator /(const RMatrix & a,double  b){
	RMatrix temp(a);
	return(temp /= b);	
}

bool RMatrix::operator ==(const RMatrix & a)const{
	if(row != a.row || col !=a.col){
		return false;
	}
	else{
		for (int i = 0; i < row; ++i) {
        	for (int j = 0; j < col; ++j) {
            	if( p[i][j] != a.p[i][j] ){
            		return false;
            	}
        	}
    	}
 		return true;
	}
}

bool RMatrix::operator !=(const RMatrix & a)const{
	return !(*this == a);
}

RMatrix RMatrix::transpose(){
	RMatrix temp(col,row);
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            temp.p[j][i] = p[i][j];
        }
    }
    return temp;
}

bool RMatrix::is_square()const{
	if (row != col){
		return false;
	}
	return true;
}

bool RMatrix::is_diagonal()const{
	if (!is_square()){
		return false;
	}
	else{
		for (int i = 0; i < row; ++i) {
        	for (int j = 0; j < col; ++j) {
            	if( (p[i][j] != 0) && (i != j) ){
            		return false;
            	}
        	}
   		}
 	return true;
	}
}

bool RMatrix::is_zero()const{
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if( p[i][j] != 0 ){
            	return false;
            }
        }
    }
 	return true;
}
bool RMatrix::is_unit()const{
	if ( (!is_square()) || (!is_diagonal()) ){
		return false;
	}
	for (int i = 0; i < row; ++i){
		if(p[i][i] != 1){
			return false;
		}
	}
	return true;
}

bool RMatrix::is_symmetric()const{
	RMatrix temp(*this);
	if(*this == temp.transpose()){
		return true;
	}
	return false;
}

bool RMatrix::is_upper_triangular()const{
	if (!is_square()){
		return false;
	}
	for (int i = 0; i < row; ++i) {
        for (int j = 0; j<i; ++j) {
            if( p[i][j] != 0 ){
            	return false;
            }
        }
    }
 	return true;

}

bool RMatrix::is_lower_triangular()const{
	RMatrix temp(*this);
	temp = temp.transpose();
	if (temp.is_upper_triangular()){
		return true;
	}
	return false;
}

double RMatrix::determ(double** a, int n){
	int i,j,j1,j2;
	double det = 0;
	double **m = NULL;

	if (n == 1){
		det = a[0][0];
	} 
   	else{
   		if (n == 2){
   			det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   		} 
   		else{
      		det = 0;
      		for (j1=0;j1<n;j1++){
         		m = new double*[n-1];
         		for (i=0;i<n-1;i++){
           			m[i] = new double[n-1];
         		}
         		for (i=1;i<n;i++){
            		j2 = 0;
            		for (j=0;j<n;j++){
               			if (j == j1)
                  		continue;
               			m[i-1][j2] = a[i][j];
               			j2++;
           			 }
         		}
         		det += pow(-1.0,1.0+j1+1.0) * a[0][j1] * determ(m,n-1);
         		for (i=0;i<n-1;i++){
           			delete[] m[i];
         		}
         		delete[] m;
      		}
   		}
   }
   return det;
}

double RMatrix::determinant(){
	if (is_square()){
		return determ(this->p,row);
	}
	else{
		throw "Can't calculate determinant";
	}
}

void RMatrix::cofactor(double** a,int n ,double** b){
	int i,j,ii,jj,i1,j1;
	double det;
	double **c;

   	c = new double*[n-1];
   	for (i=0;i<n-1;i++){
     	c[i] = new double[n-1];
   	}

   	for (j=0;j<n;j++) {
      	for (i=0;i<n;i++) {
        	i1 = 0;
         	for (ii=0;ii<n;ii++){
            	if (ii == i)
               		continue;
            	j1 = 0;
            	for (jj=0;jj<n;jj++) {
              		if (jj == j)
                  		continue;
               		c[i1][j1] = a[ii][jj];
               		j1++;
            	}
            i1++;
         	}
        det = determ(c,n-1);

        b[i][j] = pow(-1.0,i+j+2.0) * det;
      	}
  	}
   	for (i=0;i<n-1;i++){
       	delete[] c[i];
    }
	delete[] c;
}

RMatrix RMatrix::inverse(){ //пусть сама ругается на нулевой определитель
	if (is_square()){
		if(double det = determinant()){
			RMatrix B(row,row);
			cofactor(this->p,row,B.p);
			B = (1/det)* B.transpose();
			return B;			
		}
	}
	else{throw "Can't calculate inversion";}
}


RMatrix RMatrix::power(int n){
	if (is_square()){
		if (n==0){
			RMatrix temp(row,row);
			for (int i = 0; i < row; ++i){
				temp.p[i][i] = 1;
			}
			return temp;	
		}
		else{
			RMatrix temp(*this);
			if(n<0){
				temp=temp.inverse();
			}
			RMatrix temp2(temp);
			for(int j =0;j < abs(n)-1;++j){
				temp = temp * temp2;		
			}
			return temp;
		}

	}
	else{throw "Can't calculate power";}	
}

double RMatrix::norm()const{
	double m=0;
	for(int i=0;i<row;++i){
		double sum = 0;
		for(int j=0;j<col;++j){
			sum += abs(p[i][j]);
		}
		m=max(sum,m);
	}
	return m;
}






int main(){
	cout<<"Hello!"<<endl;
	for(;;){
		cout<<"What would you like to do? (Choose the number)"<<endl;
		cout<<"0 -  Matrix addition"<<endl;
		cout<<"1 -  Matrix subtraction"<<endl;
		cout<<"2 -  Matrix multiplication"<<endl;
		cout<<"3 -  Matrix multiplication by number"<<endl;
		cout<<"4 -  Matrix division by number"<<endl;
		cout<<"5 -  Matrix inversion"<<endl;
		cout<<"6 -  Matrix transposition"<<endl;
		cout<<"7 -  Matrix determinant"<<endl;
		cout<<"8 -  Matrix norm"<<endl;
		cout<<"9 -  Matrix power"<<endl;
		cout<<"10 - Matrix type: square;"<<endl<<"                  diagonal;"<<endl<<"                  zero;"<<endl<<"                  unit;"<<endl<<"                  symmetric;"<<endl<<"                  upper triangular;"<<endl<<"                  lower triangular"<<endl;
		
		int flag;
		cin>>flag;
		if(flag == 0 || flag == 1){
			int n,m;
			cout<<"Enter matrix size (n x m):"<<endl;
			cout<<"n = "; cin>>n;
			cout<<"m = "; cin>>m;
			RMatrix M1(n,m), M2(n,m);
			cout<<"Enter M1"<<endl; cin>>M1;
			cout<<"Enter M2"<<endl; cin>>M2;
			cout<<"Result is:"<<endl<<M1 + M2 * pow(-1,flag) <<endl;
		}
		if (flag == 2){
			int n1,m1,n2,m2;
			cout<<"Enter M1 matrix size (n1 x m1):"<<endl;
			cout<<"n1 = "; cin>>n1;
			cout<<"m1 = "; cin>>m1;
			cout<<"Enter M2 matrix size (n2 x m2):"<<endl;
			cout<<"n2 = "; cin>>n2;
			cout<<"m2 = "; cin>>m2;
			RMatrix M1(n1,m1), M2(n2,m2);
			cout<<"Enter M1"<<endl; cin>>M1;
			cout<<"Enter M2"<<endl; cin>>M2;
			cout<<"Result is:"<<endl<<M1 * M2 <<endl;
		}
		if (flag == 3 ||flag == 4){
			int n,m;
			double x;
			cout<<"Enter matrix size (n x m):"<<endl;
			cout<<"n = "; cin>>n;
			cout<<"m = "; cin>>m;
			RMatrix M1(n,m);
			cout<<"Enter M"<<endl; cin>>M1;
			cout<<"Enter number"<<endl; cin>>x;
			if(flag == 3){
				cout<<"Result is:"<<endl<<M1*x<<endl;
			}
			else{
				cout<<"Result is:"<<endl<<M1/x<<endl;
			}
		}
		if (flag == 5 || flag ==  7 || flag == 9){
			int n;
			cout<<"Enter matrix size (n x n):"<<endl;
			cout<<"n = "; cin>>n;
			RMatrix M1(n,n);
			cout<<"Enter M"<<endl; cin>>M1;
			switch(flag){
				case 5: cout<<"Result is:"<<endl<<M1.inverse()<<endl; break;
				case 7: cout<<"Result is:"<<endl<<M1.determinant()<<endl; break;
				case 9:{int x;
						cout<<"Enter exponent"<<endl; cin>>x;
						cout<<"Result is:"<<endl<<M1.power(x)<<endl;
						}
						break;
			}
		}

		if (flag == 6 || flag ==10 || flag ==8){
			int n,m;
			cout<<"Enter matrix size (n x m):"<<endl;
			cout<<"n = "; cin>>n;
			cout<<"m = "; cin>>m;
			RMatrix M1(n,m);
			cout<<"Enter M"<<endl; cin>>M1;
			if(flag == 6){
				cout<<"Result is:"<<endl<<M1.transpose()<<endl;
			}
			if(flag ==10){
				cout<<"Result is:"<<endl;
				if (M1.is_square()){cout<< "is square"<<endl;} else {cout<<"not square"<<endl;}
				if (M1.is_diagonal()){cout<< "is diagonal"<<endl;}else {cout<<"not diagonal"<<endl;}
				if (M1.is_symmetric()){cout<< "is symmetric"<<endl;}else {cout<<"not symmetric"<<endl;}
				if (M1.is_unit()){cout<< "is unit"<<endl;}else {cout<<"not unit"<<endl;}
				if (M1.is_zero()){cout<< "is zero"<<endl;}else {cout<<"not zero"<<endl;}
				if (M1.is_upper_triangular()){cout<< "is upper triangular"<<endl;}else {cout<<"not upper triangular"<<endl;}
				if (M1.is_lower_triangular()){cout<< "is lower triangular"<<endl;}else {cout<<"not lower triangular"<<endl;}
			}
			if(flag == 8){
				cout<<"Result is:"<<endl<<M1.norm()<<endl;
			}

		}

		cout<<"Anything else?..(yes/no)"<<endl;
		string s;
		cin >> s;

		if (s=="no"){cout<<"See you later (◕ ◡ ◕) "<<endl; break;};

	}
	return 0;
}