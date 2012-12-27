/** Simple code to test the BossMask class
 *  and make sure it works as specified.
 *
 *  Nikhil Padmanabhan
 */

#include <iostream>
#include "BossMask.h"

using namespace std;
using namespace Boss;

int main() {
	BossMask mask("bossmask_template.cfg");
	mask.print();
}
