// Struct para um vector bidimensional; declarar com vec2<tipo>
// Usar para vectores (matemáticos) em vez de usar vectors (do C++)
// Mais rápido e com sintaxe melhor (v.x,v.y em vez de v[0],v[1])
template <class T>
struct vec2 {
    T x,y;
};

#include <fstream>
#include <cmath>
#include <cstdlib> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <stack>
#include "srj.hpp"

#define alfa 1

using namespace std;

// Atalho para vector 2d de doubles: vec2d
typedef vec2<double> vec2d;

// Usar como atalho para x*x, equivalente
template <typename T> inline T sq(T x) {
    return x*x;
}

/*
const int    Lx = 128, Ly = 128;  // Dimensoes da grelha
const double medini = -0.2, ge = 0.126;
const double Amp = 10.55;  // Amplitude da forca
const double D = 11;
const double dt = 0.02;  // Passo de tempo
const int    rad = 5;  // Raio das celulas
const double vbase = 0.01;
const double valoralfa = 0.065;
const int    nmax = 4;  // Numero de tip cells maximo
const double rho0,L0,M,vconc;
*/

int p_j=0, dp_n;
double dp_res;

vector<double> srj;
vector<vec2<double>> tips;

//double vmax,Pmax;
double a[Lx][Ly],w[Lx][Ly],csi[Lx][Ly],v[Lx][Ly],a_med;

double p_res[Lx][Ly];
int p_n[Lx][Ly];

// Condições de fronteira periódicas para x e y
inline int bx(int xx) {
    return (xx+Lx)%Lx;
}
inline int by(int yy) {
    return (yy+Ly)%Ly;
}

int t;

int main()
{
    int passo;
    double rand;

    void ini();
    void step();
    void out(double x, double y);
    void outint(int ff);
    void csicalc(double raio, int index);
    void newtip();

    double find(double inic, double cy);
    double medp(vec2<double>);

    int neigh(int i,int j);

    vec2<double> findxy(vec2<double> pos, vec2<double> gradxy);
    vec2<double> gradxy(vec2<double> V); 

    int count_chunks();

    ofstream resultados("res");

    random_device seed;
    mt19937_64 rand_gen(seed());
    uniform_real_distribution<double> dist01(0.0,1.0);

    passo = 500;
    const int passotips = 500;
    int nchunks = 1;

    srj = schedule<5>(get_scheme<5>(Lx));

	for (int iN=0;iN<iNf;iN++) {
		//vmax=0.3;
		//Pmax=0.03;

		ini();

		for (t=0;t<tf;t++) {
			step();
			if ((t+100) % passotips	== 0) {
				printf("		t = %d\n",t);
				printf("         x        y       ∇x       ∇y\n");
				printf("*****************************************\n");
				for (int k=0;k<tips.size();k++) {
					const auto grad = gradxy(tips[k]);
					tips[k] = findxy(tips[k],grad);
					if (tips[k].x > Lx-1-rad || tips[k].x < rad) {
						printf("Tip %d out of bounds.\n",k+1);
						return -1;
					}
					printf("Tip %d %-8.4lf %-8.4lf %-8.4lf %-8.4lf\n",k+1,tips[k].x,tips[k].y,grad.x,grad.y);
				}
				printf("*****************************************\n");
				printf('\n');

				if (count_chunks() > nchunks) {
					printf("Vasos partidos!\n\n");
					nchunks++;
				}

				csicalc(double(rad), 0);
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
				printf("(Novo output)\n\n");
			if ((t+2) % passo == 0)
				outint(t+2);
		}
		resultados << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
		cout       << Pmax << "\t" << vmax << "\t" << a_med << "\t" << t*dt << "\n";
	}
}

// Retorna o número de vasos distintos
int count_chunks()
{
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
			stack<vec2<int>> s;
			s.push(vec2<int> {i,j});
			while (!s.empty()) {
				auto top = s.top();
				s.pop();
				if (!V[top.x+1][top.y] && a[top.x+1][top.y] > 0) {
					s.push(vec2<int> {top.x+1,top.y});
					V[top.x+1][top.y] = 1;
				}
				if (!V[top.x-1][top.y] && a[top.x-1][top.y] > 0) {
					s.push(vec2<int> {top.x-1,top.y});
					V[top.x-1][top.y] = 1;
				}
				if (!V[top.x][top.y+1] && a[top.x][top.y+1] > 0) {
					s.push(vec2<int> {top.x,top.y+1});
					V[top.x][top.y+1] = 1;
				}
				if (!V[top.x][top.y-1] && a[top.x][top.y-1] > 0) {
					s.push(vec2<int> {top.x,top.y-1});
					V[top.x][top.y-1] = 1;
				}
			}
		}
	}
	return r;
}

// Retorna a média pesada
double medp(vec2<double> V)
{
    const int x_up = ceil(V.x);
    const int x_down = floor(V.x);

    const int y_right = ceil(V.y);
    const int y_left = floor(V.y);

    const double dxup=x_up-V.x;
    const double dyr=y_right-V.y;
    const double dxdown=x_down-V.x;
    const double dyl=y_left-V.y;

    const double dist1=dxup*dxup+dyl*dyl;
    const double dist2=dxup*dxup+dyr*dyr;
    const double dist3=dxdown*dxdown+dyl*dyl;
    const double dist4=dxdown*dxdown+dyr*dyr;

    const double average = (a[x_up][y_left]/dist1+a[x_up][y_right]/dist2+a[x_down][y_left]/dist3+a[x_down][y_right]/dist4);
    const double sum = 1/dist1+1/dist2+1/dist3+1/dist4;

    return average/sum;
}


// Retorna o índice do elemento com valor máximo num vector
int maxval(const vector<double>& list)
{
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

bool neigh(int i, int j)
{
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

// Calcula a posição da nova tip e adiciona ao vector das tips
void newtip()
{
	int maxval(const vector<double>& list);
	double phimed;
	vector<vec2<int>> interface ((Lx-2)*(Ly-2), {0,0});
	vector<double> vegf ((Lx-2)*(Ly-2), 0.);
	vec2<int> pick, point;

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

	int senti = 1;

	while (senti!=0) {
		int maxvegf = maxval(vegf);
		pick = interface[maxvegf];

		for (int i=0;i<tips.size();i++) {
			const double dist1 = sqrt((pick.x-tips[i].x)*(pick.x-tips[i].x)+(pick.y-tips[i].y)*(pick.y-tips[i].y));
			const double dist2 = sqrt((pick.x-tips[i].x-Lx)*(pick.x-tips[i].x-Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    
			const double dist3 = sqrt((pick.x-tips[i].x+Lx)*(pick.x-tips[i].x+Lx)+(pick.y-tips[i].y)*(pick.y-tips[i].y));    

			if (dist1<(4.*rad) || dist2<(4.*rad) || dist3<(4.*rad)) {
				vegf[maxvegf]=-1000;
				break;
			}
			if (i==tips.size()-1)
				senti = 0;
		}
	}
	const vec2<double> ret {double(pick.x),double(pick.y)};
	tips.push_back(ret);
}

// Imprime os outputs para os ficheiros
void outint (int ff)
{
	char s[32];
	ofstream wout;
	ofstream aout;
	ofstream vout;
	//ofstream csiout;
	ofstream tout;
	double prolif(int i, int j);
	void prolifupdate();

	sprintf(s,"data/wout_%d",ff);
	wout.open(s);
	sprintf(s,"data/aout_%d",ff);
	aout.open(s);
	sprintf(s,"data/vout_%d",ff);
	vout.open(s);
	//sprintf(s,"data/csiout_%d",ff);
	//csiout.open(s);
	sprintf(s,"data/tout_%d",ff);
	tout.open(s);

	for (int j=0;j<Ly;j++){
		for (int i=0;i<Lx;i++){
			wout<<w[i][j]<<" ";
			aout<<a[i][j]<<" ";
			vout<<v[i][j]<<" ";
			//csiout<<csi[i][j]<<" ";
		}
		wout<<"\n";
		aout<<"\n";
		vout<<"\n";
		//csiout<<"\n";
	}
	for (int i=0; i<tips.size(); i++) {
		tout << tips[i].x << " " << tips[i].y << "\n";
	}

	aout.close();
	wout.close();
	vout.close();
	//csiout.close();
	tout.close();
}

// Inicializa os arrays e as tips
void ini()
{
	a_med=0;
	for (int i=0;i<Lx;i++){
		for (int j=0;j<Ly;j++){
			a[i][j] = i>10 && i<60 ? 1 : -1;
			w[i][j] = 0;
			csi[i][j] = 0;
			i>60?(v[i][j] = vconc):(v[i][j]==0);
			a_med += a[i][j];
		}
	}
	a_med /= (Lx*Ly);

	ifstream tips_in ("tips.in");
	double x,y;
	while (tips_in >> x >> y) {
		tips.push_back({x,y});
	}
}

// Efectua um passo
void step()
{
    void poisson();
    void ch();

    poisson();
    ch();
}

inline double f(double z) {
    return -valoralfa*z;
}

void ch()
{
	double Q[Lx][Ly],mu[Lx][Ly],an[Lx][Ly],vn[Lx][Ly],aloc;
	double I1[Lx][Ly],I2[Lx][Ly],I3[Lx][Ly],IE[Lx][Ly];
	double prolif(int i, int j);
	void prolifupdate();
	double consumo(int i, int j);

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

	for(int i=0;i<Lx;i++) {
		for(int j=0;j<Ly;j++) {
			aloc=a[i][j];
			mu[i][j]=rho0*(-aloc+aloc*aloc*aloc-(a[bx(i+1)][j]+a[bx(i-1)][j]+a[i][by(j-1)]+a[i][by(j+1)]-4*aloc))-ge*Q[i][j]+valoralfa*csi[i][j]/L0;
		}
	}

	a_med=0;
	prolifupdate();
	for(int i=0;i<Lx;i++) {
		for(int j=0;j<Ly;j++) {
			an[i][j]=a[i][j]+dt*(M1*(mu[bx(i+1)][j]+mu[bx(i-1)][j]+mu[i][by(j-1)]+mu[i][by(j+1)]-4*mu[i][j])+(t*dt>10?prolif(i,j):0) );
#if (alfa !=0)
			an[i][j]=an[i][j]+2*ge*IE[i][j]*valoralfa*dt/L0;
#endif
			a_med+=an[i][j];
		}
	}
	a_med/=(Lx*Ly);


	for(int i=1;i<Lx-1;i++) {
		for(int j=0;j<Ly;j++) {
			vn[i][j]=v[i][j]+D*dt*(v[i+1][j]+v[i-1][j]+v[i][by(j+1)]+v[i][by(j-1)]-4*v[i][j]-consumo(i,j));
		}
	}

	for(int j=0;j<Ly;j++){
		vn[0][j]=vn[1][j];
		vn[Lx-1][j]=vn[Lx-2][j];
	}

	for(int i=0;i<Lx;i++){
		for(int j=0;j<Ly;j++){
			a[i][j]=an[i][j];
			v[i][j]=vn[i][j];
		}
	}

}

void poisson()
{
	double wn[Lx][Ly];
	const double tol = 1E-4;
	const double dtau0 = 0.24;

	double diff=1;
	do {
		for (int i=0; i<srj.size() && diff>tol; i++) {
			const double dtau = srj[i]/4;
			double sum=0;
			for (int i=0;i<Lx;i++) {
				for (int j=0;j<Ly;j++) {
					wn[i][j]=w[i][j]+dtau*(w[bx(i+1)][j]+w[bx(i-1)][j]+w[i][by(j-1)]+w[i][by(j+1)]-4*w[i][j]-(f(a[i][j])+csi[i][j])/L0);
					sum+=wn[i][j];
				}
			}
			sum/=(Lx*Ly);
			diff = 0;
			for (int i=0;i<Lx;i++) {
				for (int j=0;j<Ly;j++) {
					wn[i][j]=wn[i][j]-sum;
					diff+=fabs(wn[i][j]-w[i][j]);
					w[i][j]=wn[i][j];
				}
			}
			diff/=(Lx*Ly);
		}
	} while(diff>tol);
}

vec2<double> gradxy(vec2<double> V)
{
    vec2<double> grad {0,0};

    const int cx = lround(V.x);
    const int cy = lround(V.y);

    grad.x = (v[cx+1][cy]-v[cx-1][cy])/2.;
    grad.y = (v[cx][cy+1]-v[cx][cy-1])/2.;

    return grad;
}

void csicalc(double raio, int index)
{
	double sum;
	double Force[Lx][Ly];
	double csin[Lx][Ly];
	ofstream csiout;
	char s[20];

	for (int i=0;i<Lx;i++) {
		for (int j=0;j<Ly;j++){
			Force[i][j] = 0.;
			for (int k=0;k<tips.size();k++) {	

				const vec2<double> vgrad = gradxy(tips[k]);
				const double cos_ = vgrad.x/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));
				const double sin_ = vgrad.y/(sqrt(vgrad.x*vgrad.x+vgrad.y*vgrad.y));

				const double cx=i-tips[k].x;
				const double cy=j-tips[k].y;

				Force[i][j] += 2*Amp/raio/raio*exp(-(cx*cx+cy*cy)/raio/raio)*((1-2*cx*cx/raio/raio)*cos_+2*cx*cy/raio/raio*sin_);
			}
			csi[i][j]=csin[i][j]=0;
		}
	}

	const double tol = 1E-6;
	const double dtau0 = 0.24;

	double diff;
	do {
		for (int i=0; i<srj.size(); i++) {
			const double dtau = srj[i]/4;

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
			if (diff<tol)
				goto exit;
		}
	} while(diff>tol);
exit:

	sprintf(s,"cout.%d",index);
	csiout.open(s);

	for (int i=0;i<Lx;i++){
		for (int j=0;j<Ly;j++){
			csiout << csi[i][j] << " ";
		}
		csiout << "\n";
	}
	csiout.close();
}

vec2<double> findxy(vec2<double> pos, vec2<double> gradxy)
{
	const double modulo=sqrt(gradxy.x*gradxy.x+gradxy.y*gradxy.y);

	gradxy.x/=modulo;
	gradxy.y/=modulo;

	pos.x-=2*gradxy.x;
	pos.y-=2*gradxy.y;

	double delta=2.;
	vec2<double> posn={0,0};

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

// Calculo da proliferação com DP
double prolif(int i, int j)
{
	if (a[i][j]<=0.5) {
		return 0;
	}

	double res=0;

	int n=0;
	
	int xr = 0;
	int yr = rad;
	for (int x=0; x<















	if (j == p_j+1) {  // DP
		for (int y=j-rad; y<=j+rad; y++) {
			for (int x=i-rad;x<=i+rad;x++) {
				if (sq(x-i) + sq(y-j) <= sq(rad)) {
					res += p_res[bx(x)][by(y)];
					n += p_n[bx(x)][by(y)];
				}
			}
		}





		n = dp_n - p_n[bx(i)][by(p_j-rad)] + p_n[bx(i)][by(j+rad)];
		res = dp_res - p_res[bx(i)][by(p_j-rad)] + p_res[bx(i)][by(j+rad)];
	}
	else {
		for (int x=i-rad;x<=i+rad;x++) {
			//for (int y=j-sqrt(sq(x-i)-sq(rad);y<=j+sqrt(sq(x-i)-sq(rad));y++) {
			for (int y=j-rad; y<=j+rad; y++) {
				if (sq(x-i) + sq(y-j) <= sq(rad)) {
					res += p_res[bx(x)][by(y)];
					n += p_n[bx(x)][by(y)];
				}
			}
		}
	}
	dp_n = n;
	dp_res = res;
	p_j = j;
	return n>0 ? res/n : 0;
}

void prolifupdate()
{
	for (int i=0; i<Lx; i++) {
		for (int j=0; j<Ly; j++) {

			const double ra = v[i][j];
			const double aloc = a[i][j];
			const double rc = csi[i][j]-valoralfa*aloc;

			if (aloc>0.5 && rc>-valoralfa+0.05) {
				p_res[i][j] = ra>vmax ? Pmax : ra*Pmax/vmax;
				p_n[i][j] = 1;
			}
			else {
				p_res[i][j] = 0;
				p_n[i][j] = 0;
			}
		}
	}
}

double consumo(int i, int j)
{
	return (a[i][j]>0 ? (a[i][j]<1 ? 0.1*a[i][j] : 0.1) : 0)*v[i][j];
}
