#include "BCH.h"

vector<vector<bool>> prim_poly = {
	{1,1}, {1,1,0}, {1,1,0,0},{1,0,1,0,0},{1,1,0,0,0,0},{1,0,0,1,0,0,0},{1,0,1,1,1,0,0,0}
};

BCH::BCH(int n, int d)
{
	if (n < d || d < 5 || log2(n + 1) != (int)log2(n + 1))
	{
		cout << "ERROR" << endl;
		exit(0);
	}
	this->n = n;
	this->d = d;
	this->m = (int)log2(n + 1);
	this->t = (d - 1) / 2;

	cout << "Galois field by modulo 2^" << m << endl;
	this->gf = create_GF(m, prim_poly[m - 2]);
	printGF(gf, m);

	this->cyclotomic_classes = create_cyclotomic_classes(n);
	cout << "Cyclotomic classes by modulo 2^" << m << endl;
	print(cyclotomic_classes);


	this->g = create_g();
	cout << "Generating polynomial: " << endl;
	print(g);

	this->k = this->n - deg(g);

	cout << "n: " << n << endl;
	cout << "k: " << k << endl;
	cout << "d: " << d << endl;
	cout << "t: " << t << endl << endl;

}

vector<vector<bool>> BCH::create_GF(int m, vector<bool>& primitive_poly) {
	vector<vector<bool>> mtrx;
	vector<bool> tmp(m);
	mtrx.push_back(tmp); //add 0,0,..,0
	tmp[0] = 1;
	mtrx.push_back(tmp); //add 1,0,..,0
	for (int i = 2; i < pow(2, m); i++) {
		tmp.clear();
		if (mtrx[i - 1][m - 1] == 1) {
			for (int j = 0; j < m; j++) {
				if (j == 0)
					tmp.push_back((primitive_poly[j] + 0) % 2);
				if (j > 0)
					tmp.push_back((primitive_poly[j] + mtrx[i - 1][j - 1]) % 2);
			}
		}
		else {
			for (int j = 0; j < m; j++) {
				if (j == 0)
					tmp.push_back(0);
				if (j > 0)
					tmp.push_back(mtrx[i - 1][j - 1]);
			}
		}
		mtrx.push_back(tmp);
	}
	return mtrx;
}

vector<vector<int>> BCH::create_cyclotomic_classes(int mod)
{
	vector<vector<int>> mtrx;
	vector<bool> pin(mod);
	vector<int> tmp;
	int conjugate = 0;
	bool flag = true;
	while (flag)
	{
		tmp.push_back(conjugate);
		conjugate = (conjugate * 2) % mod;
		pin[conjugate] = 1;
		if (conjugate == tmp[0])
		{
			mtrx.push_back(tmp);
			tmp.clear();
			int j = 0;
			while (j < mod && pin[j] == 1)
				j++;
			if (j == mod)
				flag = false;
			else conjugate = j;
		}
	}
	return mtrx;
}

vector<bool> BCH::create_g()
{
	//polynomial root selection (cyclotomic classes selection)
	vector<vector<int>> polynomial_alpha;
	int cntr = 1;
	if (d % 2 == 0 && d <= (n / 2) + 1)
		polynomial_alpha.push_back(cyclotomic_classes[0]);
	for (size_t i = 1; i < d - 1; i++)
	{
		for (size_t j = 0; j < cyclotomic_classes.size(); j++)
		{
			if (find(cyclotomic_classes[j].begin(), cyclotomic_classes[j].end(), cntr) != cyclotomic_classes[j].end() &&
				find(polynomial_alpha.begin(), polynomial_alpha.end(), cyclotomic_classes[j]) == polynomial_alpha.end())
			{
				polynomial_alpha.push_back(cyclotomic_classes[j]);
				break;
			}

		}
		cntr++;
	}

	cout << "Polynomial roots: " << endl;
	print(polynomial_alpha);

	// finding primitive polynomials by cyclotomic roots
	// f = (x - a^i1)(x - a^i2)... i beloning cyclotomic class
	vector<vector<bool>> vectors_sub_g;
	for (size_t i = 0; i < polynomial_alpha.size(); i++)
		vectors_sub_g.push_back(get_min_poly(polynomial_alpha[i]));

	//sub_polynomial multiplication, result = polynom G
	vector<bool> res = vectors_sub_g[0];
	for (size_t i = 1; i < vectors_sub_g.size(); i++)
		res = multiplication(res, vectors_sub_g[i]);

	return res;
}

vector<bool> BCH::get_min_poly(vector<int> cyclotomic_class)
{
	vector<vector<bool>> multiplier(2);
	multiplier[0] = gf[1];
	vector<vector<bool>> res;
	res.push_back(gf[1]);
	res.push_back(gf[cyclotomic_class[0] + 1]);
	for (size_t i = 1; i < cyclotomic_class.size(); i++)
	{
		multiplier[1] = gf[cyclotomic_class[i] + 1];
		res = multiplication(res, multiplier);
	}
	vector<bool> itogo;
	for (size_t i = 0; i < res.size(); i++)
		itogo.push_back(res[i][0]);
	
	return itogo;
}

vector<vector<bool>> BCH::multiplication(vector<vector<bool>> a, vector<vector<bool>> b)
{
	vector<vector<bool>> res(a.size() + b.size() - 1);
	for (size_t i = 0; i < a.size(); i++)
		for (size_t j = 0; j < b.size(); j++)
			res[i + j] = sub_xor(res[i + j], sub_multiply(a[i], b[j]));
	return res;
}

vector<bool> BCH::multiplication(vector<bool> a, vector<bool> b)
{
	vector<bool> res(a.size() + b.size() - 1);
	for (size_t i = 0; i < a.size(); i++)
		for (size_t j = 0; j < b.size(); j++)
			res[i + j] = (res[i + j] + a[i] * b[j]) % 2;
	return res;
}

vector<bool> BCH::sub_multiply(vector<bool> a, vector<bool> b)
{
	vector<vector<bool>>::iterator it = find(gf.begin(), gf.end(), a);
	int index = it - gf.begin() - 1;//power alpha
	it = find(gf.begin(), gf.end(), b);
	index += it - gf.begin() - 1;//summ powers
	index %= n;
	return gf[index + 1];
}

vector<bool> BCH::sub_xor(vector<bool> a, vector<bool> b)
{
	if (a.size() == 0)
		a.resize(this->m);
	for (int i = 0; i < a.size(); i++)
		a[i] = a[i] ^ b[i];
	return a;
}

/**
 * Modulo division
 * @param v - binary vector (polynomial)
 * @param h - binary vector (polynomial)
 * @return binary vector (polynomial) - remainder of the division
 */
vector<bool> BCH::division(vector<bool> v, vector<bool> h)
{
	vector<bool> g;
	if (v.size() < h.size() && vector_weight(v) < vector_weight(h))
		return v;

	// need vectors of the same size
	// for example: v = 1101; h = 01 -> vector h becomes 0001
	for (int i = 0; i < v.size() - h.size(); i++)
		g.push_back(0);

	for (int i = 0; i < h.size(); i++)
		g.push_back(h[i]);

	while (true)
	{
		int i, j = 0;
		for (i = 0; i < v.size(); i++)
			if (v[i] == 1) break;

		for (j = 0; j < g.size(); j++)
			if (g[j] == 1) break;

		if (i > j)
			break;

		for (int k = j; k < g.size(); k++)
		{
			v[i] = v[i] ^ g[j];
			i++;
			j++;
		}
		i = 0;
		j = 0;
	}
	return v;
}

/**
 * @param v - binary vector
 * @return vector weight - for the case of binary vectors weight is equal to the number of units in the vector
 */
int BCH::vector_weight(vector<bool> v)
{
	int res = 0;
	for (int i = 0; i < v.size(); i++)
		if (v[i] == 1) res++;

	return res;
}

/**
 * @param v - binary vector (polynomial)
 * @return polynomial degree
 */
int BCH::deg(vector<bool> v)
{
	int size_r = 0;
	for (int i = 0; i < v.size(); i++)
		if (v[i] == 1) {
			size_r = v.size() - i - 1;
			break;
		}
	return size_r;
}

vector<bool> BCH::sindrom(vector<bool> b, vector<bool> g)
{
	return division(b, g);
}

vector<bool> BCH::coding(vector<bool>m, vector<bool> g)
{
	return multiplication(m, g);
}

vector<bool> BCH::getG()
{
	return g;
}



void print(vector<bool> v)
{
	for (size_t j = 0; j < v.size(); j++)
		cout << v[j] << " ";
	cout << endl << endl;
}

void print(vector<vector<int>> v)
{
	for (size_t i = 0; i < v.size(); i++) {
		for (size_t j = 0; j < v[i].size(); j++)
			cout << v[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void print(vector<vector<bool>> v)
{
	for (size_t i = 0; i < v.size(); i++) {
		for (size_t j = 0; j < v[i].size(); j++)
			cout << v[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void printGF(vector<vector<bool>> gf, int m)
{
	for (int i = 0; i < pow(2, m); i++) {
		if (i <= 1)
			cout << i << "\t";
		else cout << "a^" << i - 1 << "\t";
		for (int j = 0; j < m; j++)
			cout << gf[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}