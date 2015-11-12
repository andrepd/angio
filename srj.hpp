#include <vector>

using namespace std;

template <const int P>
struct scheme {
	//const int p = P;
	const int N;
	const double O[P];
	const int q[P];
	//int M;
	//double OM[P];
};

/*
vector<scheme<2>> schemes = {
	{128,{425.8,0.9742},{1,130}},
};
*/

vector<scheme<5>> schemes = {
	//{64,{1228.8,220.14,26.168,3.1668,0.63890}
	{128,{4522.0,580.86,50.729,4.5018,0.67161},{1,3,16,73,250}},//,343,{343,114,21,4}}
	{256,{16459.,1513.4,97.832,6.4111,.70531},{1,4,26,142,605}},
	{512,{59226.0,3900.56,187.53,9.1194,0.73905},{1,6,40,277,1500}},
};

template<const int P>
scheme<P> get_scheme(int L) {
	if (L <= schemes.front().N)
		return schemes.front();
	for (int i=1; i<schemes.size(); i++) {
		if (schemes[i-1].N < L && L <= schemes[i].N)
			return schemes[i];
	}
	return schemes.back();
}

template<const int P>
vector<double> schedule(const scheme<P>& s) {
	int M = 0;
	for (int i=0; i<P; i++)
		M += s.q[i];
	vector<double> r (M,0);
	for (int i=0; i<P; i++) {
		double pos = i+8;
		for (int j=0; j<s.q[i]; j++) {
			//cerr << i << " " << j << " " << (int)pos << "\n";
			while (r[(int)pos%M] != 0)
				pos += 1.;
			r[(int)pos%M] = s.O[i];
			pos += (double)M/s.q[i];
		}
	}

	// debug
	/*
	vector<int> count (P,0);
	for (auto i: r) {
		for (int j=0; j<P; j++) {
			if (i == s.O[j]) {
				count[j]++;
				break;
			}
		}
	}
	
	for (auto i: count)
		cerr << i << " ";
	cerr << "\n";
	*/

	return r;
}