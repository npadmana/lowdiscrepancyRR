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
	BossMask mask("boss_survey.ply",0.0, true);
	mask.print();
}
