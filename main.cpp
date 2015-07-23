#include <fstream>
//#include <math.h>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib> 
//#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
//#include <omp.h>
#include <random>
//#include <chrono>
#include <stack>

#define alfa 1
#define teste 1

using namespace std;

// Struct para um vector bidimensional; declarar com vect2<tipo>
// Usar para vectores (matemáticos) em vez de usar vectors (do C++)
// Mais rápido e com sintaxe melhor (v.x,v.y em vez de v[0],v[1])
template <class T>
struct vect2 {
	T x,y;
};

inline double sq(double x) {
	return x*x;
}

const int    Lx = 128, Ly = 128;  // Dimensoes da grelha
const double medini = -0.2, ge = 0.126;
const double Amp = 10.55;  // Amplitude da forca
const double D = 11;
const double dt = 0.02;  // Passo de tempo
const int    rad = 5;  // Raio das celulas
const double vbase = 0.01;
const double valoralfa = 0.065;
const int    nmax = 4;  // Numero de tip cells maximo

//vector<vector<double>> fontes={/*{120,60},{120,70},{120,80},{120,90}*/};//fontes de VEGF

vector<vect2<double>> tips;  // Lista das tip cell

double vmax,Pmax;
double a[Lx][Ly],w[Lx][Ly],csi[Lx][Ly],v[Lx][Ly],a_med;

int bx(int xx){
	return (xx+Lx)%Lx;
}

int by(int yy){
	return (yy+Ly)%Ly;
}

int t;
int main() {
	double vmenor,vmaior,Pmenor,Pmaior,dv,dP;
	int passo;
	double rand;

	void ini();
	void step();
	void out(double x, double y);
	void outint(int ff);
	void csicalc(const vector<vect2<double>>& tips, double raio, int index);
	void newtip();

	double find(double inic, double cy);
	double medp(vect2<double>);

	int neigh(int i,int j);

	vect2<double> findxy(vect2<double> pos, vect2<double> gradxy);
	vect2<double> gradxy(vect2<double> V); 
	//vect2<double> grad {0,0};

	int count_chunks();

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
	uniform_real_distribution<double> dist01(0.0,1.0);

	//srand48(time(0));
	passo = 500;
	const int passotips = 500;
	const int iNf=1, tf=100000;
	int nchunks = 1;
	
	for (int iN=0;iN<iNf;iN++) {
		vmax=0.3;
		Pmax=0.03;

		ini();

		for (t=0;t<tf;t++) {
			//cerr << "t=" << t;
			step();
			//cerr << "  ";
			//for (auto i: tips) {
			//	cerr << i.x << "," << i.y << " ";
			//}
			//cerr << "\n";
			//cerr << " OK\n";
			if ((t+100) % passotips	== 0) {
				cout << "t=" << t << '\n';
				for (int k=0;k<tips.size();k++) {
					const auto grad = gradxy(tips[k]);
					tips[k] = findxy(tips[k],grad);
					//tips[k][0]=find(tips[k][0]-2,tips[k][1]);
					if (tips[k].x > Lx-1-rad || tips[k].x < rad) {
						cout << "Tip " << k+1 << " out of bounds.\n";
						return -1;
					}
					cout << "Tip " << k+1 << " " 
						<< tips[k].x << " " << tips[k].y << " "
						<< grad.x << " " << grad.y << '\n';
				}
				cout<<'\n';

				//cout << count_chunks() << "\n";
				if (count_chunks() > nchunks) {
					cout << "Vasos partidos!\n\n";
					nchunks++;
				}

				csicalc(tips, double(rad), 0);
				if (fabs(a_med)>20) 
					break;
			}
			if ((t+100) % passotips == 0 && tips.size()<nmax) {
				//rand=drand48();
				//rand = dist01(rand_gen);
				//if(rand>0.5) 
					newtip();
			}
			if ((t+1) % passo == 0)
				cout << "(Novo output)\n\n";
#if(teste==1)
			if ((t+2) % passo == 0)
				outint(t+2);
#endif
		}
		resultados << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
		cout       << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
	}
}

int count_chunks() {
	bool V[Lx][Ly];
	for (int i=0; i<Lx; i++) {
		for (int j=0; j<Ly; j++) {
			V[i][j] = 0;
		}
	}
	int r = 0;
	for (int i=0; i<Lx; i++) {
		for (int j=0; j<Ly; j++) {
			if (V[i][j] || a[i][j] < 0)
				continue;
			r++;
			stack<vect2<int>> s;
			s.push(vect2<int> {i,j});
			while (!s.empty()) {
				auto top = s.top();
				s.pop();
				if (!V[top.x+1][top.y] && a[top.x+1][top.y] > 0) {
					s.push(vect2<int> {top.x+1,top.y});
					V[top.x+1][top.y] = 1;
				}
				if (!V[top.x-1][top.y] && a[top.x-1][top.y] > 0) {
					s.push(vect2<int> {top.x-1,top.y});
					V[top.x-1][top.y] = 1;
				}
				if (!V[top.x][top.y+1] && a[top.x][top.y+1] > 0) {
					s.push(vect2<int> {top.x,top.y+1});
					V[top.x][top.y+1] = 1;
				}
				if (!V[top.x][top.y-1] && a[top.x][top.y-1] > 0) {
					s.push(vect2<int> {top.x,top.y-1});
					V[top.x][top.y-1] = 1;
				}
			}
		}
	}
	return r;
}

double medp(vect2<double> V) {
	double sum = 0;
	double sumw = 0;

	const int cx = lround(V.x);
	const int cy = lround(V.y);

	//for (int i=xc-1;i<(xc+rad+1);i++) {
	//	for (int j=yc-1;j<(yc+rad+1);j++) {
	for (int i=cx-rad;i<=cx+rad;i++) {
		for (int j=cy-rad;j<=cy+rad;j++) {
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
	int res = 0;
	double maxval = 0;

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

	// TODO
	int senti = 1;

	while (senti!=0) {
		int maxvegf = maxval(vegf);
		pick = interface[maxvegf];

		for (int i=0;i<tips.size();i++) {
			const double dist1 = sqrt((pick.x-tips[i].x)*(pick.x-tips[i].x)+(pick.y-tips[i].y)*(pick.y-tips[i].y));
			const double dist2 = sqrt((pick.x-tips[i].x-Lx)*(pick.x-tips[i].x-Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    
			const double dist3 = sqrt((pick.x-tips[i].x+Lx)*(pick.x-tips[i].x+Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    

			if (dist1<(4.*rad) || dist2<(4.*rad) || dist3<(4.*rad)) {
				vegf[maxvegf]=-1000;  // TODO
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
	char s[32];
	ofstream wout;
	ofstream aout;
	ofstream vout;
	ofstream csiout;
	ofstream tout;
	double prolif(int i, int j);

	sprintf(s,"data/wout_%d",ff);
	wout.open(s);
	sprintf(s,"data/aout_%d",ff);
	aout.open(s);
	sprintf(s,"data/vout_%d",ff);
	vout.open(s);
	sprintf(s,"data/csiout_%d",ff);
	csiout.open(s);
	sprintf(s,"data/tout_%d",ff);
	tout.open(s);

	for (int j=0;j<Ly;j++){
		for (int i=0;i<Lx;i++){
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
	for (int i=0; i<tips.size(); i++) {
		tout << tips[i].x << " " << tips[i].y << "\n";
	}
	aout.close();
	wout.close();
	vout.close();
	csiout.close();
	tout.close();
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
	a_med=0;
	for (int i=0;i<Lx;i++){
		for (int j=0;j<Ly;j++){
			a[i][j] = i>Lx/5 && i<Lx/5+10 ? 1 : -1;
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

inline double f(double z){
	return -valoralfa*z;
}

void ch() {
	//double f(double z);
	double Q[Lx][Ly],mu[Lx][Ly],an[Lx][Ly],vn[Lx][Ly],aloc;
	double I1[Lx][Ly],I2[Lx][Ly],I3[Lx][Ly],IE[Lx][Ly];
	double prolif(int i, int j);
	double consumo(int i, int j);

	//double maiorv,maiora,maiorw;

	for (int i=0;i<Lx;i++) {
		for (int j=0;j<Ly;j++) {
			aloc=a[i][j];
			const double txx=w[bx(i+1)][j]+w[bx(i-1)][j]-2*w[i][j];
			const double tyy=w[i][by(j+1)]+w[i][by(j-1)]-2*w[i][j];
			const double txy=(w[bx(i+1)][by(j+1)]-w[bx(i-1)][by(j+1)]+w[bx(i-1)][by(j-1)]-w[bx(i+1)][by(j-1)])/4.0;
			Q[i][j]=txx*txx+tyy*tyy+2*txy*txy-(f(aloc)+csi[i][j])*(f(aloc)+csi[i][j])/2.0;
#if (alfa != 0)
			I1[i][j]=aloc*(txx-(f(aloc)+csi[i][j])/2.0);
			I2[i][j]=aloc*(tyy-(f(aloc)+csi[i][j])/2.0);
			I3[i][j]=aloc*txy;
#endif
		}
	}

#if (alfa != 0)
	for (int i=0;i<Lx;i++) {
		for (int j=0;j<Ly;j++) {
			IE[i][j]=(I1[bx(i+1)][j]+I1[bx(i-1)][j]-2*I1[i][j])+(I2[i][by(j+1)]+I2[i][by(j-1)]-2*I2[i][j])+2*(I3[bx(i+1)][by(j+1)]-I3[bx(i-1)][by(j+1)]+I3[bx(i-1)][by(j-1)]-I3[bx(i+1)][by(j-1)])/4.0;
		}
	}
#endif

	for(int i=0;i<Lx;i++){
		for(int j=0;j<Ly;j++){
			aloc=a[i][j];
			mu[i][j]=-aloc+aloc*aloc*aloc-(a[bx(i+1)][j]+a[bx(i-1)][j]+a[i][by(j-1)]+a[i][by(j+1)]-4*aloc)-ge*Q[i][j]+valoralfa*csi[i][j];
		}
	}

	a_med=0;
	for(int i=0;i<Lx;i++){
		for(int j=0;j<Ly;j++){
			an[i][j]=a[i][j]+dt*(mu[bx(i+1)][j]+mu[bx(i-1)][j]+mu[i][by(j-1)]+mu[i][by(j+1)]-4*mu[i][j]+ (t*dt>10?prolif(i,j):0) );
#if (alfa !=0)
			an[i][j]=an[i][j]+2*ge*IE[i][j]*valoralfa*dt;
#endif
			a_med+=an[i][j];
		}
	}
	a_med/=(Lx*Ly);


	for(int i=1;i<Lx-1;i++){
		for(int j=0;j<Ly;j++){
			vn[i][j]=v[i][j]+D*dt*(v[i+1][j]+v[i-1][j]+v[i][by(j+1)]+v[i][by(j-1)]-4*v[i][j]-consumo(i,j));
		}
	}

	for(int j=0;j<Ly;j++){
		vn[0][j]=0;
		vn[Lx-1][j]=1;
	}

	for(int i=0;i<Lx;i++){
		for(int j=0;j<Ly;j++){
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
	//double f(double z);
	double wn[Lx][Ly];
	const double tol = 1E-4;
	const double dtau = 0.24;

	double diff;
	//double diff_;
	int qq = 0;
	do {
		qq++;
		double sum=0;
		for (int i=0;i<Lx;i++) {
			for (int j=0;j<Ly;j++) {
				wn[i][j]=w[i][j]+dtau*(w[bx(i+1)][j]+w[bx(i-1)][j]+w[i][by(j-1)]+w[i][by(j+1)]-4*w[i][j]-(f(a[i][j])+csi[i][j]));
				sum+=wn[i][j];
			}
		}
		sum/=(Lx*Ly);
		//diff_ = diff;
		diff = 0;
		for (int i=0;i<Lx;i++) {
			for (int j=0;j<Ly;j++) {
				wn[i][j]=wn[i][j]-sum;
				diff+=fabs(wn[i][j]-w[i][j]);
				w[i][j]=wn[i][j];
			}
		}
		diff/=(Lx*Ly);
		//int foo;
		//if (fabs(diff-diff_) > 1) cin >> foo; 
		//cerr << "    " << diff << "\n";
	} while(diff>tol);
	cerr << ">P " << qq << "\n";
	//if (qq > 16) {
		//int foo;
		//cin >> foo;
	//}
}

vect2<double> gradxy(vect2<double> V) {
	vect2<double> grad {0,0};

	const int cx = lround(V.x);
	const int cy = lround(V.y);

	grad.x = (v[cx+1][cy]-v[cx-1][cy])/2.;
	grad.y = (v[cx][cy+1]-v[cx][cy-1])/2.;

	return grad;
}

void csicalc(const vector<vect2<double>>& tips, double raio, int index) {
	double sum;
	double Force[Lx][Ly];
	//double Force_[Lx][Ly];
	double csin[Lx][Ly];
	ofstream csiout;
	char s[20];
	//double cx=0,cy=0,cynew,cxnew;
	//vect2<double> vgrad {0,0};
	//double dir=0;

	for (int i=0;i<Lx;i++) {
		for (int j=0;j<Ly;j++){
			Force[i][j] = 0.;
			//Force_[i][j] = 0.;
			for (int k=0;k<tips.size();k++) {	

				const vect2<double> vgrad = gradxy(tips[k]);
				//dir=atan(vgrad.y/vgrad.x);
				const double cos_ = vgrad.x/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));
				const double sin_ = vgrad.y/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));

				//if ( (vgrad.x>0 && vgrad.y<0) || (vgrad.x<0 && vgrad.y<0) ) 
				//	dir+=M_PI;

				const double cx=i-tips[k].x;
				const double cy=j-tips[k].y;

				//cxnew=cx*cos(dir)-cy*sin(dir);
				//cynew=cx*sin(dir)+cy*cos(dir);

				//Force_[i][j]+=-Amp*exp(-cynew*cynew/raio/raio)*exp(-cxnew*cxnew/raio/raio)*(4*cxnew*cxnew-2*raio*raio)/raio/raio/raio/raio;	
				//Force[i][j] += -Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((2+4*cx*cx/raio/raio)*cos(dir)+4*cx*cy/raio/raio*sin(dir));
				//Force[i][j] += -Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((2+4*cx*cx/raio/raio)*cos_+4*cx*cy/raio/raio*sin_);
				Force[i][j] += 2*Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((1-2*cx*cx/raio/raio)*cos_+2*cx*cy/raio/raio*sin_);
				//cerr << i << " " << j << " " << Force[i][j] << " " << Force_[i][j] << "\n";
			}
			csi[i][j]=csin[i][j]=0;
		}
	}

	const double tol=1E-6;
	const double dtau=0.24;

	double diff;
	int qq = 0;
	do {
		qq++;
		// TODO
		//diff = 0;
		for(int i=1;i<Lx-1;i++){
			for(int j=1;j<Ly-1;j++){
				csin[i][j]=csi[i][j]+dtau*(csi[bx(i+1)][j]+csi[bx(i-1)][j]+csi[i][by(j-1)]+csi[i][by(j+1)]-4*csi[i][j]-Force[i][j]);
				diff+=fabs(csin[i][j]-csi[i][j]);
			}
		}
		for(int i=0;i<Lx;i++){
			for(int j=0;j<Ly;j++){
				csi[i][j]=csin[i][j];
			}
		}
		diff/=(Lx*Ly);
		//cerr << "  " << diff << "\n";
	} while(diff>tol);
	cerr << ">C " << qq << "\n";
	//if (qq > 16) {
	//	int foo;
	//	cin >> foo;
	//}

	sprintf(s,"cout.%d",index);
	csiout.open(s);

	for (int i=0;i<Lx;i++){
		for (int j=0;j<Ly;j++){
			csiout << csi[i][j] << " ";
		}
		csiout << "\n";
	}
	csiout.close();
	//cerr << "DONE\n";
}

vect2<double> findxy(vect2<double> pos, vect2<double> gradxy) {
	const double modulo=sqrt(gradxy.x*gradxy.x+gradxy.y*gradxy.y);

	gradxy.x/=modulo;
	gradxy.y/=modulo;

	pos.x-=2*gradxy.x;
	pos.y-=2*gradxy.y;

	double delta=2.;
	vect2<double> posn={0,0};

	while (delta>0.001) {
		posn.x=pos.x+delta*gradxy.x;
		posn.y=pos.y+delta*gradxy.y;

		if (medp(posn)>0) {
			pos=posn;
		} else {
			delta/=2.;
		}
	}
	return pos;
}

// TODO
/*
double find(double inic, double cy){
	int res1,res2,med;

	res1=(int)inic;
	res2=(int)Lx;
	while(abs(res1-res2)>2){
		med=(int)(res1+res2)/2.0+1;
		if(a[res1][(int)cy]*a[med][(int)cy]<0)res2=med;
		else res1=med;
	}

	const double resultado=res1-(res2-res1)*a[res1][(int)cy]/(a[res2][(int)cy]-a[res1][(int)cy]);

	return resultado;
}
*/

double prolif(int i, int j){
	if (a[i][j]<=0.5) 
		return 0;
	else {
		double res=0;
		int n=0;
		for (int x=i-rad;x<=i+rad;x++) {
			for (int y=j-rad;y<=j+rad;y++) {
				if (sq(x-i)+sq(y-j) <= sq(rad)) {
					const double ra = v[bx(x)][by(y)];
					const double aloc = a[bx(x)][by(y)];
					const double rc = csi[bx(x)][by(y)]-valoralfa*aloc;
					if (aloc>0.5 && rc>-valoralfa+0.05) {
						res += ra>vmax ? Pmax : ra*Pmax/vmax;
						n++;
					}
				}
			}
		}
		return n>0 ? res/n : 0;
	}
}

double consumo(int i, int j) {
	return (a[i][j]>0 ? (a[i][j]<1 ? 0.1*a[i][j] : 0.1) : 0)*v[i][j];
}
