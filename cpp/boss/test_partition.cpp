/** Develop and test a code to partition a number
 *  for parallel work
 *
 *  Nikhil Padmanabhan, Yale, Dec 28, 2012
 */

#include <iostream>
#include <cmath>
#include <tuple>

using namespace std;

tuple<int, int> partition(int n) {
	// Short circuit the case of n=1
	if (n==1) return make_tuple(1,1);

	int sn = static_cast<int>(floor(sqrt(n)+0.01)); // shift slightly to avoid missing perfect squares.
	while ((n % sn)!=0) sn--;

	return make_tuple(sn, n/sn);
}

int main() {
	int n, n1, n2;
	cout << "Enter an integer :";
	cin >> n;
	tie(n1, n2) = partition(n);

	cout << n1 << " " << n2 << endl;
}
