#include "BCH.h"



int main() 
{
	BCH code(15,8);
	vector<bool> mes = {1,0,1};
	//vector<bool> a = { 1,1,0,0,1,0,0,0,1,1,1,1,0,1,0 };
	vector<bool> a = code.coding(mes, code.getG());
	cout << "Vector a:" << endl;
	print(a);
	vector<bool> e = {1,0,0,0,1,0,0,0 };
	vector<bool> b = vector_b(a, e);
	cout << "Vector b = a + e:" << endl;
	print(b);
	vector<bool> res = code.decoding(b);
	cout << "Decoded received vector:" << endl;
	print(res);
	cout << "Vector a:" << endl;
	print(a);

	return 0;
}

