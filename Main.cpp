#include "BCH.h"

vector<bool> vector_b(vector<bool> a, vector<bool> e);

int main() 
{
	BCH code(15,7);
	vector<bool> a = { 1,1,0,0,1,0,0,0,1,1,1,1,0,1,0 };
	//vector<bool> a = code.coding(mes, code.getG());
	cout << "Vector a:" << endl;
	print(a);
	vector<bool> e = { 0,0,0,0,1,0,0,0,0,0,0,1,0,0,0 };
	vector<bool> b = vector_b(a, e);
	cout << "Vector b = a + e:" << endl;
	print(b);


	return 0;
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