#include "BCH.h"

vector<vector<bool>> prim_poly = {
	{1,1}, {1,1,0}, {1,1,0,0},{1,0,1,0,0},{1,1,0,0,0,0},{1,0,0,1,0,0,0},{1,0,1,1,1,0,0,0}
};

/* public foo */

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

	this->cyclotomic_classes_by_field = create_cyclotomic_classes(n);
	cout << "Cyclotomic classes by modulo 2^" << m << endl;
	print(cyclotomic_classes_by_field);


	this->g = create_g();
	cout << "Generating polynomial: " << endl;
	print(g);

	this->k = this->n - deg(g);

	cout << "n: " << n << endl;
	cout << "k: " << k << endl;
	cout << "d: " << d << endl;
	cout << "t: " << t << endl << endl;

}

vector<bool> BCH::coding(vector<bool>m, vector<bool> g)
{
	//cout << deg(g) << " - degree of the generating polynomial (deg(g) = r)" << endl;
	vector<bool> c = check_sum(m, g);
	//for (int l = 0; l < c.size(); l++)
	//	cout << c[l] << ' ';
	//cout << " - check sum" << endl;

	vector<bool> a = vector_a(m, c);
	//for (int l = 0; l < a.size(); l++)
	//	cout << a[l] << ' ';
	//cout << " - vector a" << endl;
	return a;
	//return multiplication(m, g);
}

vector<bool> BCH::getG()
{
	return g;
}

vector<bool> BCH::decoding(vector<bool> &vec)
{
	vector<vector<bool>> sndrm_cmpnnts = get_components_of_syndrome(vec);
	cout << "Components of syndrome:" << endl;
	print(sndrm_cmpnnts);
	int cntr = 0;
	for (size_t i = 0; i < sndrm_cmpnnts.size(); i++)
		if (sndrm_cmpnnts[i] != gf[0])
			cntr++;
	if (cntr == 0)
	{
		cout << "Recieved message accepted without error " << endl << endl;
		return vec;
	}

	int v = t;
	vector<vector<vector<bool>>> M;
	vector<vector<bool>>tmp;
	vector<bool> det;
	while (v > 0) {
		M.resize(v);
		tmp.resize(v);
		for (size_t i = 0; i < v; i++)
			M[i] = tmp;

		for (size_t i = 0; i < v; i++)
			for (size_t j = 0; j < v; j++)
				M[i][j] = sndrm_cmpnnts[i + j];
		cout << "Matrix component of the syndrome:" << endl;
		print(M);

		det = determinant(M,v);
		cout << "Determinant of matrix size " << v << ": ";
		if (det == gf[0]) cout << "0" << endl << endl;
		else
		{
			vector<vector<bool>>::iterator it = find(gf.begin(), gf.end(), det);
			int index = it - gf.begin() - 1;//power alpha
			cout << "alpha^" << index << endl << endl;
		}
		if (det != gf[0]) break;
		v--;
	}
	cout << "Number of errors in the received word >= " << v << endl << endl;
	cout << "Matrix M:" << endl;
	print(M);

	// Inverse Matrix (by algebraic addition)
	vector<vector<vector<bool>>> res(M.size());
	vector<vector<bool>> sub_res(M.size());
	vector<vector<vector<bool>>> temp(M.size() - 1);
	vector<vector<bool>> sub_temp(M.size() - 1);
	for (size_t i = 0; i < temp.size(); i++)
		temp[i] = sub_temp;
	for (size_t i = 0; i < res.size(); i++)
		res[i] = sub_res;
	if (M.size() == 1)
		res[0][0] = inverse_element(M[0][0], det);
	else {
		for (size_t i = 0; i < M.size(); i++) {
			for (size_t j = 0; j < M.size(); j++) {
				int n = 0;
				int m = 0;
				for (size_t k = 0; k < temp.size(); k++) {
					for (size_t l = 0; l < temp.size(); l++) {
						if (n == i)
							n++;
						if (m == j)
							m++;
						temp[k][l] = M[n][m];
						m++;
					}
					n++;
					m = 0;
				}
				res[j][i] = determinant(temp, temp.size());
			}
		}
		vector<bool> num = inverse_element(det, det);
		for (int i = 0; i < M.size(); i++) {
			for (int j = 0; j < M.size(); j++) {
				//res[i][j] = (num * res[i][j]) % module;
				res[i][j] = sub_multiply(num, res[i][j]);
			}
		}
	}
	M = res;

	// Inverse Matrix + inverse element(det(M))
	for (size_t i = 0; i < M.size(); i++)
		for (size_t j = 0; j < M[0].size(); j++)
			M[i][j] = inverse_element(M[i][j],det);
	cout << "Matrix M^-1:" << endl;
	print(M);

	// Finding the coefficient matrix of a polynomial
	vector<vector<bool>> vector_s;
	for (size_t i = 0; i < v; i++)
		vector_s.push_back(sndrm_cmpnnts[v + i]);
	vector<vector<bool>> coef_vector(M.size());
	vector<bool> coef_alpha;
	for (size_t i = 0; i < M.size(); i++)
	{
		coef_vector[i] = gf[0];
		for (size_t j = 0; j < M[0].size(); j++)
			//coef_vector[i] = coef_vector[i] + M[i][j] * vector_s[j];
			coef_vector[i] = sub_xor(coef_vector[i], sub_multiply(M[i][j], vector_s[j]));
	}
	//print(coef_vector);
	vector<vector<bool>> res_invers;
	vector<bool> is_res_invers = gf[1];
	// Chen' procedure
	vector<vector<bool>>::iterator it;
	int index;
	for (size_t i = 1; i < gf.size(); i++)
	{
		is_res_invers = gf[1];
		it = find(gf.begin(), gf.end(), gf[i]);
		for (size_t j = 0; j < coef_vector.size(); j++)
		{
			index = it - gf.begin() - 1; //power alpha
			index = (index * (coef_vector.size() - j)) % n;
			is_res_invers = sub_xor(is_res_invers, sub_multiply(coef_vector[j], gf[index + 1]));
		}
		if (gf[0] == is_res_invers)
			res_invers.push_back(gf[i]);
	}
	if (res_invers.size() < M.size())
	{
		cout << endl;
		cout << "Decoding denial, number of errors in the received word > "<< t << endl;
		exit(0);
	}

	// Getting error locators
	vector<vector<bool>>error_locators;
	for (size_t i = 0; i < res_invers.size(); i++)
		error_locators.push_back(inverse_element(gf[1], res_invers[i]));
	//print(error_locators);

	// Error polynom
	vector<int> result;
	for (size_t i = 0; i < error_locators.size(); i++)
	{
		it = find(gf.begin(), gf.end(), error_locators[i]);
		index = it - gf.begin() - 1; //power alpha
		result.push_back(index);
	}
	cout << "Vector e(x) = ";
	sort(result.begin(), result.end());
	reverse(result.begin(), result.end());
	for (size_t i = 0; i < result.size(); i++)
	{
		cout << "x^" << result[i];
		if (i != result.size() - 1)
			cout << " + ";
		else cout << endl;
	}
	cout << endl;

	// Result
	reverse(vec.begin(), vec.end());
	for (size_t i = 0; i < result.size(); i++)
	{
		if (vec[result[i]] == 1)
			vec[result[i]] = 0;
		else if (vec[result[i]] == 0)
			vec[result[i]] = 1;
	}
	reverse(vec.begin(), vec.end());

	return vec;
}

/* private foo */

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

	int cntr = 1;
	if (d % 2 == 0 && d <= (n / 2) + 1)
		polynomial_roots.push_back(cyclotomic_classes_by_field[0]);
	for (size_t i = 1; i < d - 1; i++)
	{
		for (size_t j = 0; j < cyclotomic_classes_by_field.size(); j++)
		{
			if (find(cyclotomic_classes_by_field[j].begin(), cyclotomic_classes_by_field[j].end(), cntr) != cyclotomic_classes_by_field[j].end() &&
				find(polynomial_roots.begin(), polynomial_roots.end(), cyclotomic_classes_by_field[j]) == polynomial_roots.end())
			{
				polynomial_roots.push_back(cyclotomic_classes_by_field[j]);
				break;
			}

		}
		cntr++;
	}

	cout << "Polynomial roots: " << endl;
	print(polynomial_roots);

	// finding primitive polynomials by cyclotomic roots
	// f = (x - a^i1)(x - a^i2)... i beloning cyclotomic class
	vector<vector<bool>> vectors_sub_g;
	for (size_t i = 0; i < polynomial_roots.size(); i++)
		vectors_sub_g.push_back(get_min_poly(polynomial_roots[i]));

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
	if (a == gf[0] || b == gf[0])
	{
		return gf[0];
	}
	vector<vector<bool>>::iterator it = find(gf.begin(), gf.end(), a);
	int index = it - gf.begin() - 1; //power alpha
	it = find(gf.begin(), gf.end(), b);
	index += it - gf.begin() - 1; //summ powers
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

vector<bool> BCH::inverse_element(vector<bool> chisl, vector<bool> znamen)
{
	vector<vector<bool>>::iterator it = find(gf.begin(), gf.end(), chisl);
	int index = it - gf.begin() - 1;//power alpha
	it = find(gf.begin(), gf.end(), znamen);
	index -= it - gf.begin() - 1;//summ powers
	if (index < 0)
		index += n;

	return gf[index + 1];
}

vector<vector<bool>> BCH::get_components_of_syndrome(vector<bool> b)
{
	vector<vector<int>> sub_res;
	vector<int>tmp;
	// i - alphas from 1 to 2t	
	for (size_t i = 1; i <= 2 * t; i++)
	{
		tmp.clear();
		for (size_t j = 0; j < b.size(); j++) 
		{
			if (b[j] == 1) 
				tmp.push_back(((b.size() - 1 - j) * i) % n);
		}
		sub_res.push_back(tmp);
	}
	cout << "In-between calculations components of syndrome:" << endl;
	print(sub_res);

	vector<vector<bool>> res;
	for (size_t i = 0; i < sub_res.size(); i++)
	{
		vector<bool> tmp_bool = gf[sub_res[i][0] + 1];
		for (size_t j = 1; j < sub_res[i].size(); j++)
			tmp_bool = sub_xor(tmp_bool, gf[sub_res[i][j] + 1]);

		res.push_back(tmp_bool);
	}

	return res;
}

/* helper foo */

vector<bool> BCH::determinant(vector<vector<vector<bool>>> &mas, int m) {
	int i, j, n;
	vector<vector<vector<bool>>>p(m);
	vector<vector<bool>>tmp(m);
	for (i = 0; i < m; i++)
		p[i] = tmp;
	j = 0; 
	vector<bool> d = gf[0];
	n = m - 1;
	if (m < 1) cout << "The determinant is impossible to calculate!";
	if (m == 1) {
		d = mas[0][0];
		return(d);
	}
	//if (m == 2) {
	//	//d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
	//	d = sub_xor(sub_multiply(mas[0][0], mas[1][1]), sub_multiply(mas[1][0], mas[0][1]));
	//	return(d);
	//}
	//if (m > 2) {
		for (i = 0; i < m; i++) {
			GetMatr(mas, p, i, 0, m);
			//d = d + mas[i][0] * Determinant(p, n);
			d = sub_xor(d, sub_multiply(mas[i][0], determinant(p, n)));

		}
	//}
	return(d);
}

/**
 * Getting the matrix without the i-th row and j-th column
 * @param mas - 
 * @param p -
 * @param i -
 * @param j -
 * @param m -
 */
void BCH::GetMatr(vector<vector<vector<bool>>> &mas, vector<vector<vector<bool>>> &p, int i, int j, int m) {
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki < m - 1; ki++) { // row index check
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) { // column index check
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}

/**
 * Modulo division
 * @param v - binary vector (polynomial)
 * @param h - binary vector (polynomial)
 * @return binary vector (polynomial) - remainder of the division
 */
vector<bool> division(vector<bool> v, vector<bool> h)
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
int vector_weight(vector<bool> v)
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
int deg(vector<bool> v)
{
	int size_r = 0;
	for (int i = 0; i < v.size(); i++)
		if (v[i] == 1) {
			size_r = v.size() - i - 1;
			break;
		}
	return size_r;
}

vector<bool> sindrom(vector<bool> b, vector<bool> g)
{
	return division(b, g);
}

vector<bool> vector_a(vector<bool> m, vector<bool> c)
{
	for (int i = 0; i < c.size(); i++)
		m.push_back(c[i]);

	return m;
}

vector<bool> vector_b(vector<bool> a, vector<bool> e)
{
	if (e.size() > a.size())
		throw runtime_error("vector_size(e) > vector_size(a) \n");
	while (e.size() != a.size())
		e.insert(e.begin(), 0);
	for (int i = 0; i < a.size(); i++)
		a[i] = a[i] ^ e[i];
	return a;
}

vector<bool> check_sum(vector<bool> m, vector<bool> g)
{
	vector<bool> c = m;
	vector<bool> res;

	for (int i = 0; i < deg(g); i++)
		c.push_back(0);

	c = division(c, g);

	// the size of the vector c must be equal to deg(g) 
	// for example if result of division 01 but deg(g) = 4, then next code make it to 0001)
	for (int i = 0; i < deg(g); i++)
		res.push_back(c[c.size() - i - 1]);

	reverse(res.begin(), res.end());

	return res;
}

/* output foo */

void print(vector<bool> v)
{
	for (size_t j = 0; j < v.size(); j++)
		cout << v[j] << " ";
	cout << endl << endl;
}

void print(vector<int> v)
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

void print(vector<vector<vector<bool>>> v)
{
	for (size_t i = 0; i < v.size(); i++)
	{
		for (size_t j = 0; j < v[i].size(); j++)
		{
			for (size_t k = 0; k < v[i][j].size(); k++)
			{
				cout << v[i][j][k];
			}
			cout << " ";
		}
		cout << endl;
	}
	//cout << endl;
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
