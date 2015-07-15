#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <cstdlib> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <random>
#include <chrono>

#define alfa 1
#define teste 1

using namespace std;

template <class T>
struct vect2 {
	T x,y;
};

const int Lx=128,Ly=128;  // Dimensoes da grelha
const double medini=-0.2,ge=0.126;
const double Amp=10.55;  // Amplitude da forca
const double D=11;
const double dt=0.02;  // Passo de tempo
const int rad=5;  // Raio das celulas
const double vbase=0.01;
const double valoralfa=0.065;

int nmax=2;  // Numero de tip cells maximo

//vector<vector<double>> fontes={/*{120,60},{120,70},{120,80},{120,90}*/};//fontes de VEGF

vector<vect2<double>> tips;  // Lista das tip cell

double vmax,Pmax;
double a[Lx][Ly],w[Lx][Ly],csi[Lx][Ly],v[Lx][Ly],a_med;

int t;
int bx(int xx);
int by(int yy);

int bx(int xx){
	return ((xx+Lx)%Lx);
}

int by(int yy){
	return ((yy+Ly)%Ly);
}

int main() {
	double vmenor,vmaior,Pmenor,Pmaior,dv,dP;
	int passo,iN;
	double rand;

	void ini();
	void step();
	void out(double x,double y);
	void outint(int ff);
	void csicalc(const vector<vect2<double>>& tips, double raio, int index);
	void newtip();

	double find(double inic,double cy);
	double medp(vect2<double>);

	int neigh(int i,int j);

	vect2<double> findxy(vect2<double> pos, vect2<double> gradxy);
	vect2<double> gradxy(vect2<double> V); 
	vect2<double> grad {0,0};

	ofstream resultados("res");

	tips.push_back({Lx/5+10,Ly/2.});

#if (teste == 0)
	vmenor=0.20;vmaior=1.0;
	Pmenor=0.002;Pmaior=0.01;
	dv=0.1;dP=0.001;
	passo=5000;
#endif

#if(teste == 1)
	vmenor=0.1;vmaior=0.1;
	Pmenor=0.02;Pmaior=0.02;
	dv=0.1;dP=0.01;
	passo=1000;
#endif

#if (teste == 2)
	vmenor=0.20;vmaior=1.0;
	Pmenor=0.002;Pmaior=0.01;
	dv=(vmaior-vmenor);
	dP=(Pmaior-Pmenor);
	passo=5000;
#endif

	random_device seed;
	mt19937_64 rand_gen(seed());
	uniform_real_distribution<double> dist(0.0,1.0);

	//srand48(time(0));
	passo=1000;

	for (iN=0;iN<4;iN++) {
		vmax=0.3;
		Pmax=0.03;

		ini();

		for (t=0;t<100000;t++) {
			//cerr << "t=" << t;
			step();
			//cerr << " OK\n";
			if ((t+100) % 500 == 0) {
				cout << "t=" << t << '\n';
				for (int k=0;k<tips.size();k++) {
					grad = gradxy(tips[k]);
					tips[k] = findxy(tips[k],grad);
					//tips[k][0]=find(tips[k][0]-2,tips[k][1]);
					cout<<"Tip "<<k+1<<'\t'<<tips[k].x<<'\t'<<tips[k].x<<'\t'<< grad.x << '\t' << grad.y << '\n';
				}
				cout<<'\n';
				csicalc(tips,double(rad),0);
				if(fabs(a_med)>20)break;
			}
			if ((t+100) % 500 == 0 && tips.size()<nmax) {
				//rand=drand48();
				rand = dist(rand_gen);
				if(rand>0.5) newtip();
			}
			if ((t+1) % passo == 0)
				cout << "Novo output!\n";
#if(teste==1)
			if ((t+2) % passo == 0)
				outint(t+2);
#endif
		}
		resultados << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
		cout       << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
	}
}

double medp(vect2<double> V) {
	int xc = (int)V.x;
	int yc = (int)V.y;

	double sum = 0;
	double sumw = 0;

	//for (int i=xc-1;i<(xc+rad+1);i++) {
	//	for (int j=yc-1;j<(yc+rad+1);j++) {
	for (int i=xc-rad;i<=xc+rad;i++) {
		for (int j=yc-rad;j<=yc+rad;j++) {
			const double dist = sqrt((V.x-i)*(V.x-i)+(V.y-j)*(V.y-j));
			if (dist<rad) {
				sum += a[i][j]/sqrt(dist);
				sumw += 1/sqrt(dist);
			}
		}
	}
	return sum/sumw;

	//double average=a[x_up][y_left]+(x-x_up)*(a[x_down][y_left]-a[x_up][y_left])+(y-y_left)*(a[x_up][y_right]-a[x_up][y_left]);

	//double average=(a[xc+1][yc]+a[xc-1][yc]+a[xc][yc+1]+a[xc][yc-1]
	//	+0.5*a[xc+1][yc+1]+0.5*a[xc-1][yc-1]+0.5*a[xc+1][yc-1]+0.5*a[xc-1][yc+1])/6.;

	// double average=(a[xc+1][yc]+a[xc-1][yc]+a[xc][yc+1]+a[xc][yc-1])/4.;

	//return average;
}

int maxval(const vector<double>& list) {
	int res =0;
	double maxval =0;

	for (int i=0;i<list.size();i++) {
		if (list[i]>maxval) {
			maxval=list[i];
			res=i;
		}
	}

	return res;
}

// TODO
bool neigh(int i, int j) {
	int res = 0;

	if (i==0) {
		if (a[i+1][j]*a[i][j]<0) res+=1;
		else if (a[Lx-1][j]*a[i][j]<0) res+=1;
		else if (a[i][j+1]*a[i][j]<0) res+=1;
		else if (a[i][j-1]*a[i][j]<0) res+=1;
	}

	else if (i==Lx) {
		if (a[i+1][j]*a[i][j]<0) res+=1;
		else if (a[Lx-1][j]*a[i][j]<0) res+=1;
		else if (a[i][j+1]*a[i][j]<0) res+=1;
		else if (a[i][j-1]*a[i][j]<0) res+=1;
	}

	else {
		if (a[i+1][j]*a[i][j]<0) res+=1;
		else if (a[i-1][j]*a[i][j]<0) res+=1;
		else if (a[i][j+1]*a[i][j]<0) res+=1;
		else if (a[i][j-1]*a[i][j]<0) res+=1;
	}

	if (res==0 || res==4)
	   	return 0;
	else 
		return 1;
}

void newtip() {
	int maxval(const vector<double>& list);
	double phimed;
	//vector<vector<double>> interface ((Lx-2)*(Ly-2), vector<double>(2,0.));
	//vector<double> vegf ((Lx-2)*(Ly-2), 0.), pick, point;
	vector<vect2<int>> interface ((Lx-2)*(Ly-2), {0,0});
	vector<double> vegf ((Lx-2)*(Ly-2), 0.);
   	vect2<int> pick, point;

	int pos = 0;
	for (int i=1;i<Lx-1;i++) {
		for (int j=1;j<Ly-1;j++) {
			if (neigh(i,j)==1) {
				interface[pos] = {i,j};
				vegf[pos] = v[i][j];
				pos++;
			}
		}
	}

	double dist1,dist2,dist3;
	// TODO
	int senti = 1, maxvegf;

	while (senti!=0) {
		maxvegf = maxval(vegf);
		pick = interface[maxvegf];

		for (int i=0;i<tips.size();i++) {
			dist1=sqrt((pick.x-tips[i].x)*(pick.x-tips[i].x)+(pick.y-tips[i].y)*(pick.y-tips[i].y));
			dist2=sqrt((pick.x-tips[i].x-Lx)*(pick.x-tips[i].x-Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    
			dist3=sqrt((pick.x-tips[i].x+Lx)*(pick.x-tips[i].x+Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    

			if (dist1<(4.*rad) || dist2<(4.*rad) || dist3<(4.*rad)) {
				vegf[maxvegf]=-1000;
				break;
			}
			if (i==tips.size()-1)
				senti = 0;
		}
	}
	//interface.clear();
	const vect2<double> ret {double(pick.x),double(pick.y)};
	tips.push_back(ret);
}

void outint (int ff) {
	int i,j;
	char s[20];
	ofstream wout;
	ofstream aout;
	ofstream vout;
	ofstream csiout;
	double prolif(int i, int j);

	sprintf(s,"wout_%d",ff);
	wout.open(s);
	sprintf(s,"aout_%d",ff);
	aout.open(s);
	sprintf(s,"vout_%d",ff);
	vout.open(s);
	sprintf(s,"csiout_%d",ff);
	csiout.open(s);

	for (i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			wout<<w[i][j]<<" ";
			aout<<a[i][j]<<" ";
			vout<<v[i][j]<<" ";
			csiout<<csi[i][j]<<" ";
		}
		wout<<"\n";
		aout<<"\n";
		vout<<"\n";
		csiout<<"\n";
	}
	aout.close();
	wout.close();
	vout.close();
	csiout.close();
}

void out(double x, double y){
	int i,j;
	char s[20];
	ofstream wout;
	ofstream aout;
	ofstream vout;

	sprintf(s,"wout_%6.4f_%5.3f",x,y);
	wout.open(s);
	sprintf(s,"aout_%6.4f_%5.3f",x,y);
	aout.open(s);
	sprintf(s,"vout_%6.4f_%5.3f",x,y);
	vout.open(s);

	for (i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			wout<<w[i][j]<<" ";
			aout<<a[i][j]<<" ";
			vout<<v[i][j]<<" ";
		}
		wout<<"\n";
		aout<<"\n";
		vout<<"\n";
	}
	aout.close();
	wout.close();
	vout.close();
}


void ini() {
	int i,j;

	a_med=0;

	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			a[i][j] = i<Lx/5+10 && i>Lx/5 ? 1 : -1;
			w[i][j] = 0;
			csi[i][j] = 0;
			v[i][j] = ((double)i)/(Lx-1);
			a_med += a[i][j];
		}
	}
	a_med /= (Lx*Ly);
}

void step() {
	void poisson();
	void ch();

	poisson();
	ch();
}

void ch() {
	double f(double z);
	int i,j;
	double Q[Lx][Ly],mu[Lx][Ly],an[Lx][Ly],vn[Lx][Ly],aloc;
	double txx,txy,tyy;
	double I1[Lx][Ly],I2[Lx][Ly],I3[Lx][Ly],IE[Lx][Ly];
	double prolif(int i, int j);
	double consumo(int i, int j);

	double maiorv,maiora,maiorw;

	for (i=0;i<Lx;i++) {
		for (j=0;j<Ly;j++) {
			aloc=a[i][j];
			txx=w[bx(i+1)][j]+w[bx(i-1)][j]-2*w[i][j];
			tyy=w[i][by(j+1)]+w[i][by(j-1)]-2*w[i][j];
			txy=(w[bx(i+1)][by(j+1)]-w[bx(i-1)][by(j+1)]+w[bx(i-1)][by(j-1)]-w[bx(i+1)][by(j-1)])/4.0;
			Q[i][j]=txx*txx+tyy*tyy+2*txy*txy-(f(aloc)+csi[i][j])*(f(aloc)+csi[i][j])/2.0;
#if (alfa != 0)
			I1[i][j]=aloc*(txx-(f(aloc)+csi[i][j])/2.0);
			I2[i][j]=aloc*(tyy-(f(aloc)+csi[i][j])/2.0);
			I3[i][j]=aloc*txy;
#endif
		}
	}

#if (alfa != 0)
	for (i=0;i<Lx;i++) {
		for (j=0;j<Ly;j++) {
			IE[i][j]=(I1[bx(i+1)][j]+I1[bx(i-1)][j]-2*I1[i][j])+(I2[i][by(j+1)]+I2[i][by(j-1)]-2*I2[i][j])+2*(I3[bx(i+1)][by(j+1)]-I3[bx(i-1)][by(j+1)]+I3[bx(i-1)][by(j-1)]-I3[bx(i+1)][by(j-1)])/4.0;
		}
	}
#endif

	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			aloc=a[i][j];
			mu[i][j]=-aloc+aloc*aloc*aloc-(a[bx(i+1)][j]+a[bx(i-1)][j]+a[i][by(j-1)]+a[i][by(j+1)]-4*aloc)-ge*Q[i][j]+valoralfa*csi[i][j];
		}
	}

	a_med=0;
	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			an[i][j]=a[i][j]+dt*(mu[bx(i+1)][j]+mu[bx(i-1)][j]+mu[i][by(j-1)]+mu[i][by(j+1)]-4*mu[i][j]+ (t*dt>10?prolif(i,j):0) );
#if (alfa !=0)
			an[i][j]=an[i][j]+2*ge*IE[i][j]*valoralfa*dt;
#endif
			a_med+=an[i][j];
		}
	}
	a_med/=(Lx*Ly);


	for(i=1;i<Lx-1;i++){
		for(j=0;j<Ly;j++){
			vn[i][j]=v[i][j]+D*dt*(v[i+1][j]+v[i-1][j]+v[i][by(j+1)]+v[i][by(j-1)]-4*v[i][j]-consumo(i,j));
		}
	}

	for(j=0;j<Ly;j++){
		vn[0][j]=0;
		vn[Lx-1][j]=1;
	}

	for(i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			a[i][j]=an[i][j];
			v[i][j]=vn[i][j];
		}
	}

	/*for(int i=0;i<fontes.size();i++) {
	  double fontex=fontes[i][0];
	  double fontey=fontes[i][1];
	  v[(int)fontex][(int)fontey]=1.;
	  }*/
}

void poisson() {
	double f(double z);
	int i,j;
	double wn[Lx][Ly];
	double diff,dtau;
	const double tol = 1E-3;
	double sum;

	dtau=0.24;

	diff=tol+1;
	double diff_;
	while(diff>tol){
		sum=0;
		for (i=0;i<Lx;i++) {
			for (j=0;j<Ly;j++) {
				wn[i][j]=w[i][j]+dtau*(w[bx(i+1)][j]+w[bx(i-1)][j]+w[i][by(j-1)]+w[i][by(j+1)]-4*w[i][j]-(f(a[i][j])+csi[i][j]));
				sum+=wn[i][j];
			}
		}
		sum/=(Lx*Ly);
		//diff_ = diff;
		diff=0;
		for (i=0;i<Lx;i++) {
			for (j=0;j<Ly;j++) {
				wn[i][j]=wn[i][j]-sum;
				diff+=fabs(wn[i][j]-w[i][j]);
				w[i][j]=wn[i][j];
			}
		}
		diff/=(Lx*Ly);
		//int foo;
		//if (fabs(diff-diff_) > 1) cin >> foo; 
		//cerr << "    " << diff << "\n";
	}
}

vect2<double> gradxy(vect2<double> V) {

	vect2<double> grad {0,0};

	int cx=(int)V.x;
	int cy=(int)V.y;

	grad.x=(v[cx+1][cy]-v[cx-1][cy])/2.;
	grad.y=(v[cx][cy+1]-v[cx][cy-1])/2.;

	return grad;
}

void csicalc(const vector<vect2<double>>& tips, double raio, int index) {
	int i,j;
	double diff,dtau;
	double tol;
	double sum;
	double Force[Lx][Ly];
	double Force_[Lx][Ly];
	double csin[Lx][Ly];
	ofstream csiout;
	char s[20];
	double cx=0,cy=0,cynew,cxnew;
	vect2<double> vgrad {0,0};
	double dir=0;

	for(i=0;i<Lx;i++) {
		for(j=0;j<Ly;j++){
			Force[i][j] = 0.;
			Force_[i][j] = 0.;
			for(int k=0;k<tips.size();k++) {	

				vgrad=gradxy(tips[k]);
				dir=atan(vgrad.y/vgrad.x);
				const double cos_ = vgrad.x/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));
				const double sin_ = vgrad.y/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));

				if ( (vgrad.x>0 && vgrad.y<0) || (vgrad.x<0 && vgrad.y<0) ) 
					dir+=M_PI;

				cx=i-tips[k].x;
				cy=j-tips[k].y;

				cxnew=cx*cos(dir)-cy*sin(dir);
				cynew=cx*sin(dir)+cy*cos(dir);

				Force_[i][j]+=-Amp*exp(-cynew*cynew/raio/raio)*exp(-cxnew*cxnew/raio/raio)*(4*cxnew*cxnew-2*raio*raio)/raio/raio/raio/raio;	
				//Force[i][j] += -Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((2+4*cx*cx/raio/raio)*cos(dir)+4*cx*cy/raio/raio*sin(dir));
				//Force[i][j] += -Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((2+4*cx*cx/raio/raio)*cos_+4*cx*cy/raio/raio*sin_);
				Force[i][j] += -2*Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((1-2*cx*cx/raio/raio)*cos_+2*cx*cy/raio/raio*sin_);
				//cerr << i << " " << j << " " << Force[i][j] << " " << Force_[i][j] << "\n";
			}
			csi[i][j]=csin[i][j]=0;
		}
	}

	tol=1E-7;
	dtau=0.24;

	diff=tol+1;
	while(diff>tol){
		for(i=1;i<Lx-1;i++){
			for(j=1;j<Ly-1;j++){
				csin[i][j]=csi[i][j]+dtau*(csi[bx(i+1)][j]+csi[bx(i-1)][j]+csi[i][by(j-1)]+csi[i][by(j+1)]-4*csi[i][j]-Force[i][j]);
				diff+=fabs(csin[i][j]-csi[i][j]);
			}
		}
		for(i=0;i<Lx;i++){
			for(j=0;j<Ly;j++){
				csi[i][j]=csin[i][j];
			}
		}
		diff/=(Lx*Ly);
		//cerr << "  " << diff << "\n";
	}

	sprintf(s,"cout.%d",index);
	csiout.open(s);

	for (i=0;i<Lx;i++){
		for(j=0;j<Ly;j++){
			csiout<<csi[i][j]<<" ";
		}
		csiout<<"\n";
	}
	csiout.close();
	//cerr << "DONE\n";
}

vect2<double> findxy(vect2<double> pos, vect2<double> gradxy) {

	double modulo=sqrt(gradxy.x*gradxy.x+gradxy.y*gradxy.y);

	gradxy.x/=modulo;
	gradxy.y/=modulo;

	pos.x-=2*gradxy.x;
	pos.y-=2*gradxy.y;

	double delta=2.;
	vect2<double> posn={0,0};

	while (delta>0.001) {
		posn.x=pos.x+delta*gradxy.x;
		posn.y=pos.y+delta*gradxy.y;

		if(medp(posn)>0) {
			pos=posn;
		}
		else
			delta/=2.;
	}
	return pos;
}

double find(double inic, double cy){

	double resultado;
	int res1,res2,med;

	res1=(int)inic;
	res2=(int)Lx;
	while(abs(res1-res2)>2){
		med=(int)(res1+res2)/2.0+1;
		if(a[res1][(int)cy]*a[med][(int)cy]<0)res2=med;
		else res1=med;
	}

	resultado=res1-(res2-res1)*a[res1][(int)cy]/(a[res2][(int)cy]-a[res1][(int)cy]);

	return resultado;
}

double f(double z){
	return -valoralfa*z;
}

double prolif(int i, int j){

	int x,y,n;
	double res,ra,rc,aloc;

	if(a[i][j]<=0.5) res=0;
	else{
		res=0;
		n=0;
		for (x=i-rad;x<=i+rad;x++){
			for(y=j-rad;y<=j+rad;y++){
				if((x-i)*(x-i)+(y-j)*(y-j)<=rad*rad){
					ra=v[bx(x)][by(y)];
					aloc=a[bx(x)][by(y)];
					rc=csi[bx(x)][by(y)]-valoralfa*aloc;
					if(aloc>0.5 && rc>-valoralfa+0.05){
						res+=ra>vmax?Pmax:ra*Pmax/vmax;
						n++;
					}
				}
			}
		}
		res = n>0 ? res/n : 0;
	}
	return res;
}

double consumo(int i, int j){
	return (a[i][j]>0 ? (a[i][j]<1 ? 0.1*a[i][j] : 0.1) : 0)*v[i][j];
}
