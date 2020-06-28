#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class BCH {
private:
	int n, d, m, t, k;
	vector<bool>g;
	vector<vector<bool>>gf;
	vector<vector<int>> cyclotomic_classes;

public:
	BCH(int n, int d);
	vector<bool> getG();
	vector<bool> coding(vector<bool>m, vector<bool> g);
	vector<bool> decoding();

	///
	vector<vector<bool>> create_GF(int m, vector<bool>& primitive_poly);
	vector<vector<int>> create_cyclotomic_classes(int mod);
	vector<bool> create_g();
	vector<bool> get_min_poly(vector<int> cyclotomic_class);
	vector<vector<bool>> multiplication(vector<vector<bool>> a, vector<vector<bool>> b);
	vector<bool> multiplication(vector<bool> a, vector<bool> b);
	vector<bool> sub_multiply(vector<bool> a, vector<bool> b);
	vector<bool> sub_xor(vector<bool> a, vector<bool> b);

	vector<bool> sindrom(vector<bool> b, vector<bool> g);
	vector<bool> division(vector<bool> v, vector<bool> h);
	int vector_weight(vector<bool> v);
	int deg(vector<bool> v);
};
void print(vector<vector<int>> v);
void print(vector<vector<bool>> v);
void print(vector<bool> v);
void printGF(vector<vector<bool>> gf, int m);